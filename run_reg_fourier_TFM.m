% function run_reg_fourier_TFM
%RUN_REGULAR_FOURIER_TFM.M
%
% Compute tractions from substrate displacements.
% 
% This script calls reg_fourier_TFM.m written by Yunfei Huang and Benedikt 
% Sabass, and published as Sci. Rep., 2019, vol. 9, pp. 1?16. 
% https://www.nature.com/articles/s41598-018-36896-x 

% Optional: It can be useful to run batch jobs by running this script as a
% function.
% 
% This script written by Jacob Notbohm and Aashrith Saraswathibhatla
% University of Wisconsin-Madison, 2015-2020
%
% This script uses the following functions:
%  Kabsch.m
%  inpaint_nans.m
%  smooth2a.m
%  reg_fourier_TFM.m
% This script requires a file called 'ExperimentalSettings.txt.' See readme
% for more information.

clear;
close all;
clc;

%% --- USER INPUTS ---
% Name of file with displacement data
filename = 'fidic_results.mat';
% Name to save drift-corrected displacements and tractions
savename = 'tract_results.mat';
% Input name of image file containing the domain where cells are located.
% This script uses the mean displacements for pixels outside the domain to
% compute the rigid body drift; it shifts the displacements to correct for
% the drift. The domain should be an image that's the same size as the
% images used in the DIC with values of 1 inside the domain and values of
% zero outside the domain. To shift data using all displacement data,
% without excluding data points inside a domain, set domainname = [].
domainname = 'domain.tif';
% domainname = [];
% Number of time points to run. Set to empty aray [] to run all time
% points.
num_images = [];

% Optional: Crop displacement data from edges of image. This can be done if
% there's drift or other errors (e.g., due to low image quality,
% vignetting, etc.) at edges of images. To crop no data, set this parameter
% to 0.
crop_val = 10;

% Get some inputs from experimental settings file.
fid = fopen('ExperimentalSettings.txt');
txtcell = cell2mat(textscan(fid,'%f %*[^\n]')); % '%*[^\n]' skips the remainder of each line
pix_size = txtcell(1)*1e6; % Pixel size, um
E = txtcell(2);
nu = txtcell(3);
fclose(fid);


%% --- LOAD DATA ---

% Load data
load(filename);
% Number of time points
if isempty(num_images)
    % If there's a domain that's been computed, use the size of the domain
    % to get num_images
    if ~isempty(domainname)
        info = imfinfo(domainname);
        num_images = length(info);
        if num_images > size(u,3) % Check to see which is bigger: domain or number of correlations
            num_images = size(u,3);
        end
    else % Otherwise, use the number of DIC correlations
        num_images = size(u,3);
    end
end

[M, N, P] = size(u);

% Make the image square. reg_fourier_TFM.m works only for square images 
if (M>N) % if #rows are more than # colms; crop rows
    cr = round((M-N)/2);
    % Crop rows
    x = x(cr:M-cr, 1:N);
    y = y(cr:M-cr, 1:N);
    u = u(cr:M-cr, 1:N,:);
    v = v(cr:M-cr, 1:N,:); 
elseif (N>M) % if #colmns are more than #rows; crop colms
    cc = round((N-M)/2);
    % Crop columns
    x = x(1:M, cc:N-cc);
    y = y(1:M, cc:N-cc);
    u = u(1:M, cc:N-cc,:);
    v = v(1:M, cc:N-cc,:);    
end
% The above if-elseif statement will be off by one or two arrays 
% because of rounding error. Perform a final step to make the image 
% square while also cropping the image
M = min(size(x)); % number of rows before cropping
N = M;            % number of cols before cropping
% Crop data
x = x(crop_val+1:M-crop_val, crop_val+1:N-crop_val);
y = y(crop_val+1:M-crop_val, crop_val+1:N-crop_val);
u = u(crop_val+1:M-crop_val, crop_val+1:N-crop_val,:);
v = v(crop_val+1:M-crop_val, crop_val+1:N-crop_val,:);

% final size of arrays after cropping and making them square
[M, N] = size(x);

%% --- COMPUTE TRACTIONS ---

% Preallocate
u2 = zeros(M,N,P);
v2 = u2;
tx = u2;
ty = u2;
for k = 1:num_images 
    u_k = u(:,:,k);
    v_k = v(:,:,k);
    % Set nan values to zero.
    u_k(isnan(u_k))=0; v_k(isnan(v_k))=0;
    
    % Correct for image rotation/translation
    if ~isempty(domainname) % Use domain to compute image rotation and translation
        % Get k-th domain
        domain = imread(domainname,k);
        % Downsample domain
        % x and y grid points start at w0/2 and end w0/2 before the image ends.
        % First crop off edges so that domain matches start and end points of x
        % and y.
        domain = domain(min(y(:)):max(y(:)),min(x(:)):max(x(:)));
        domain = downsample(domain,d0); % downsample number of rows
        domain = downsample(domain',d0)'; % downsample number of cols
        % Set max value of domain to 1
        domain = domain/max(domain(:));
        % Get coords of locations outside domain
        out_domain = logical(1-domain);
        % Set boundaries of out_domain to zero
        out_domain(1,:) = false(1,N);
        out_domain(M,:) = false(1,N);
        out_domain(:,1) = false(M,1);
        out_domain(:,N) = false(M,1);
        % Erode out_domain using disks. Units are d0-spaced pix, so use
        % small disk radius
        SE = strel('disk',3,0);
        out_domain = imerode(out_domain,SE);
        
        % --- Correct for rigid body rotations ---
        x_rot = x(out_domain);
        y_rot = y(out_domain);
        u_rot = u_k(out_domain);
        v_rot = v_k(out_domain);
        % Arrays to compare using the Kabsch algorithm
        P = [x_rot' ; y_rot'];
        Q = [(x_rot+u_rot)' ; (y_rot+v_rot)'];
        % Get rotation using the Kabsch algorithm. U is the rotation matrix
        % about the centroid given by p0 or q0.
        [U, r, lrms, p0, q0] = Kabsch(P, Q);
        % Unrotate all data
        Q_all = [(x(:)+u_k(:))' ; (y(:)+v_k(:))'];
        v1 = ones(1,size(Q_all,2)) ;     % row vector of ones
        % Matrix U is the rotation matrix to go from P to Q. I actually
        % want to go from Q to P so switch the sign on the off diagonal
        % components.
        U(1,2) = -U(1,2);
        U(2,1) = -U(2,1);
        Q_unrotated = U*(Q_all-q0*v1) + q0*v1;
        % Displacements are rows of Q_unrotated minus the x or y
        % coordinates
        u_vector = Q_unrotated(1,:)' - x(:);
        v_vector = Q_unrotated(2,:)' - y(:);
        % Convert vector to matrix
        u_k = reshape(u_vector,size(u_k));
        v_k = reshape(v_vector,size(v_k));
        
        % --- Correct for rigid body translations ---
        % Subtract off mean displacements where out_domain==1
        u_k = u_k - mean(u_k(out_domain==1));
        v_k = v_k - mean(v_k(out_domain==1));
        
    else % Use all data
         % When location of domain isn't given, use median as an estimate 
         % of the average rigid body translation. This isn't always a good
         % idea; sometimes it will create systematic errors.
        u_k = u_k - nanmedian(u_k(:));
        v_k = v_k - nanmedian(v_k(:));
    end
    
    % If any data is large, it's likely to be an error, so interpolate its value
    umag = sqrt(u_k.^2 + v_k.^2);
    % For a normal random variable, there's a 3e-5 chance of getting
    % displacements greater than 4 standard deviations from zero, but the
    % data isn't normal, so actual probability is greater than 3e-5. Use 5
    % standard deviations on displacement magnitude as threshold.
    idx = find(abs(umag) > 10*std(umag(:)));
    % To interpolate, set to nan first; then use inpaint_nans
    u_k(idx)=nan; v_k(idx)=nan;
    u_k = inpaint_nans(u_k);
    v_k = inpaint_nans(v_k);
    
    % Filter displacements. This is optional. It will reduce spatial
    % resolution slightly, but it will also reduce noise.
    u_k = smooth2a(u_k,1); % The 1 gives a 3x3 mean smoothing
    v_k = smooth2a(v_k,1);

    % converting to microns 
    % naming this to new variable to store values in pixels
    x1 = x*pix_size;
    y1 = y*pix_size;
    u_k1 = u_k*pix_size;
    v_k1 = v_k*pix_size;

    % Mirror displacement data to top, bottom, left, and right. The final
    % size of the displacement aray will be 3x3 times larger. Mirroring
    % gets rid of high frequency noise, especially at the edges of the
    % image. It generally isn't important when signal to noise is high.
    [M,N] = size(u_k1);
    u_mirrored = zeros(3*M,3*N);
    v_mirrored = zeros(3*M,3*N);
    u_mirrored(1:M,1:N) = rot90(u_k1,2); % The '2' rotates by 90 degrees twice
    v_mirrored(1:M,1:N) = rot90(v_k1,2);
    u_mirrored(1:M,N+1:2*N) = flipud(u_k1);
    v_mirrored(1:M,N+1:2*N) = flipud(v_k1);
    u_mirrored(1:M,2*N+1:3*N) = rot90(u_k1,2);
    v_mirrored(1:M,2*N+1:3*N) = rot90(v_k1,2);
    u_mirrored(M+1:2*M,1:N) = fliplr(u_k1);
    v_mirrored(M+1:2*M,1:N) = fliplr(v_k1);
    u_mirrored(M+1:2*M,N+1:2*N) = u_k1;
    v_mirrored(M+1:2*M,N+1:2*N) = v_k1;
    u_mirrored(M+1:2*M,2*N+1:3*N) = fliplr(u_k1);
    v_mirrored(M+1:2*M,2*N+1:3*N) = fliplr(v_k1);
    u_mirrored(2*M+1:3*M,1:N) = rot90(u_k1,2);
    v_mirrored(2*M+1:3*M,1:N) = rot90(v_k1,2);
    u_mirrored(2*M+1:3*M,N+1:2*N) = flipud(u_k1);
    v_mirrored(2*M+1:3*M,N+1:2*N) = flipud(v_k1);
    u_mirrored(2*M+1:3*M,2*N+1:3*N) = rot90(u_k1,2);
    v_mirrored(2*M+1:3*M,2*N+1:3*N) = rot90(v_k1,2);
    u_k1 = u_mirrored;
    v_k1 = v_mirrored;
    
    %FFT of displacement field, which is required for reg_fourier_TFM.m
    Ftu(:,:,1) = fft2(u_k1);
    Ftu(:,:,2) = fft2(v_k1);
    
    % Combine x and y coordinates into one array called grid_mat 
	% for compatibility with reg_fourier_TFM.m
    grid_mat(:,:,1) = y1;
    grid_mat(:,:,2) = x1;
    i_max = size(u_k1,2);
    j_max = size(v_k1,1); 
    
    % Compute tractions. Inputs are um and Pa
    % inputs for reg_Fourier_TFM(fft_displ_x,fft_displ_y,reg_param,E,nu,
    %     grid_spacing (microns), #rows, #clmns, pix_size,zdepth)
    % notes: 
    % - reg_param is set to zero. Increasing this values makes calcs smoother
    % - E (modulus) is set to one. I tried changing it here but apparently
    % the equations are written for 1 value
    % - zdepth: This is to compute tractions for a plane different from the 
    % one on which displacements are measured. In our data, the fluorescent
	% particles are on the top plane so displacements are measured on the top
	% plane and zdepth=0
    [tk,~,~] = reg_fourier_TFM(Ftu(:,:,1), Ftu(:,:,2), 0, 1, nu, ...
        d0*pix_size, i_max, j_max, grid_mat, pix_size, 0); 
    %rescale Young's modulus
    tk = E*tk;
    
    % Crop mirrored tractions to convert to original size
    tk = tk(M+1:2*M,N+1:2*N,1:2);
    
    % Add data to displacement and traction arrays
    u2(:,:,k) = u_k;
    v2(:,:,k) = v_k;
    tx(:,:,k) = tk(:,:,1);
    ty(:,:,k) = tk(:,:,2);
end
% Rename arrays
u = u2;
v = v2;
% Note: Arrays u and v are corrected data (with corrections for rigid body
% motion). Original (uncorrected) arrays u and v will not be saved here.

% Save data
save(savename,'u','v','tx','ty','inc','w0','d0','x','y');


