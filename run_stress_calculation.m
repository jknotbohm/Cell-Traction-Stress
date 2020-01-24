% function run_stress_calculation(traction_data,savename,domainname,BC)
%RUN_STRESS_CALCULATION
% 
% run_stress_calculation(traction_data,savename,domainname,BC)
%
% Compute monolayer stresses. This follows the procedure the following
% papers:
%  Tambe et al, Nat Mater, 2011, 10:469-475
%  Tambe et al, Plos One, 2013, e55172
% The finite element solver is coded into a separate script called
% compute_stress.m.
%
% INPUTS
%   traction_data: name of mat file with traction data
%   savename: name of file to save with computed stresses
%   domainname: multipage tif file giving coordinates of domain on which to
%       compute stresses
%   BC: Boundary condition information. Options: 'none', 'all',
%       'left-right', 'top', 'bottom', 'left', 'right'. See documentation
%       of compute_stress.m for description
%
% OUTPUTS
%   mat file containing computed stresses with name given by savename
%
%
% Optional: It can be useful to run batch jobs by running this script as a
% function. To do this, comment out text in section USER INPUTS
%
% Written by Jacob Notbohm, University of Wisconsin-Madison, 2013-2018
%

% This script uses the following functions:
%  force_moment_equilibrium.m
%  compute_stress.mfind
% This script requires a file called 'ExperimentalSettings.txt.' See readme
% for more information.

%% --- USER INPUTS ---
% Comment these out if running as a function

clear;
close all;
clc;

% Number of timepoints to consider. To run all time points, set to []
num_images = [];
% Name of mat file with traction data
traction_data = 'displ_tractions_w0=64.mat';
% Name of file to save data
savename = 'stresses_BCleft.mat';
% Name of multipage tif file for domain
domainname = 'domain.tif';
% Boundary condition to use. Options: 'none', 'all', 'left-right', 'top',
% 'bottom', 'left', 'right'. See documentation of compute_stress.m for
% description
BC = 'left';

%% --- LOAD DATA ---

% Load data
load(traction_data);
% Total of correlations
if isempty(num_images)
    num_images = size(tx,3);
end
% Determine whether geometry is edge/strip or hole/island from
% ExperimentalSettings. txt file
fid = fopen('ExperimentalSettings.txt');
txtcell = cell2mat(textscan(fid,'%f %*[^\n]')); % '%*[^\n]' skips the remainder of each line
pix_size = txtcell(1); % m
nu = txtcell(5); % Poisson's ratio (assuming elastic solid) or dimensionless ratio of viscosities (assuming viscous fluid) for cell monolayer (-)
h = txtcell(6); % Monolayer height (m)
edge_island = txtcell(7); % Geometry of cell layer. 1 for edge/strip, 2 for island
fclose(fid);
if edge_island~=1 && edge_island~=2
    error('Line 7 of the ExperimentalSettings.txt file must be either 1 (edge/strip) or 2 (island).');
end

%% --- COMPUTE STRESSES ---

[M, N, ~] = size(tx);
Sxx = zeros(M,N,num_images);
Syy = zeros(M,N,num_images);
Sxy = zeros(M,N,num_images);

for k=1:num_images
    % Get k-th tractions
    tx_k = tx(:,:,k);
    ty_k = ty(:,:,k);
    % Get k-th domain
    domain = imread(domainname,k);
    domain = double(domain);
    domain = domain/max(domain(:));
    domain = logical(domain);
    % Downsample domain
    % x and y grid points start at w0/2 and end w0/2 before the image ends.
    % First crop off edges so that domain matches start and end points of x
    % and y.
    domain = domain(min(y(:)):max(y(:)),min(x(:)):max(x(:)));
    domain = downsample(domain,d0); % downsample number of rows
    domain = downsample(domain',d0)'; % downsample number of cols
    if any(domain(:)==1) % Only run stress computation if there are nonzero values of domain (ie, the island was found by the edge finding algorithm)
        
        % Expand domain by dilating slightly. This assures all traction
        % data underneath cells is used in the stress calculation. Too much
        % dilation is not an issue, because calcaulated stresses will be 0
        % outside the cell domain.
        SE = strel('disk',6,0);
        domain = imdilate(domain,SE);
        
        % Tractions inside domain
        tx_domain = nan*zeros(size(tx_k)); tx_domain(domain) = tx_k(domain);
        ty_domain = nan*zeros(size(ty_k)); ty_domain(domain) = ty_k(domain);
        
        % --- Prepare traction data to be input into stress calculation ---
%         % Optional: For data in edge/strip geometry, it's common for data
%         % to lie on the left side. If it lies on the right, flip it left 
%         % to right
%         if edge_island==1
%             domain_col1 = domain(:,1);
%             if mean(domain_col1) < 0.5
%                 tx_domain = fliplr(tx_domain);
%                 ty_domain = fliplr(ty_domain);
%                 domain = fliplr(domain);
%                 tx_domain = -tx_domain; % When flipping data left to right, change sign on horizontal tractions
%             end
%         end
        
        if edge_island==1 % edge/strip geometry
            % Crop off edge data to get rid of erroneous tractions at
            % boundaries
            [M, N] = size(x);
            xcrop = x(2:M-1,2:N-1);
            ycrop = y(2:M-1,2:N-1);
            tx_domain_crop = tx_domain(2:M-1,2:N-1);
            ty_domain_crop = ty_domain(2:M-1,2:N-1);
            domain_crop = domain(2:M-1,2:N-1);
            
            % Assemble data into vectors for input into FE function
            xd = xcrop(domain_crop)*pix_size; % Units: m
            yd = ycrop(domain_crop)*pix_size;
            txd = tx_domain_crop(domain_crop); % Units: Pa
            tyd = ty_domain_crop(domain_crop);
        
        elseif edge_island==2 % island geometry
            % For this geometry, we know that the island should be in force
            % and moment equilibrium. Correct for nonzero moments and
            % forces in the traction data.
            % Center of island domain
            xc = nanmean(sum(x.*domain,2)./sum(domain,2));
            yc = nanmean(sum(y.*domain,1)./sum(domain,1));
            % Get polar coordinates r and theta, where origin is center of island
            r = sqrt( (x-xc).^2 + (y-yc).^2 );
            theta = atan2(y-yc,x-xc);
            % Correct for force and moment equilibrium
            [tx_equil, ty_equil] = force_moment_equilibrium(tx_domain,ty_domain,r,theta);
            [tx_equil, ty_equil] = force_moment_equilibrium(tx_equil,ty_equil,r,theta);
            [tx_equil, ty_equil] = force_moment_equilibrium(tx_equil,ty_equil,r,theta);
            [tx_equil, ty_equil] = force_moment_equilibrium(tx_equil,ty_equil,r,theta);
            % Turn nan values into zeros for FE computation
            tx_equil(isnan(tx_equil))=0;
            ty_equil(isnan(ty_equil))=0;
            
            % Assemble data into vectors for input into FE function
            xd = x(domain)*pix_size; % Units: m
            yd = y(domain)*pix_size;
            txd = tx_equil(domain); % Units: Pa
            tyd = ty_equil(domain);
        end
        
        % Constants K1 and K2 for stress computation script are bulk and
        % shear moduli (assuming a solid) or viscosities (assuming a fluid)
        K1 = 1/(3*(1-2*nu));  % bulk
        K2 = 1/(2*(1+nu));    % shear
        % Run FE based stress computation script.
        [Sxx_k, Syy_k, Sxy_k, xc, yc, ~, ~] = compute_stress(xd,yd,txd,tyd, K1,K2, h,BC);
        
        % Interpolate stresses to a grid
        Sxx(:,:,k) = griddata(xc,yc,Sxx_k,x*pix_size,y*pix_size);
        Syy(:,:,k) = griddata(xc,yc,Syy_k,x*pix_size,y*pix_size);
        Sxy(:,:,k) = griddata(xc,yc,Sxy_k,x*pix_size,y*pix_size);
        
        disp(['Stress computation for time ',num2str(k),' of ',num2str(num_images),' complete.'])
    else
        disp(['Unable to compute stresses for time ',num2str(t)'.'])
    end
    
end

% Set nan values to 0
Sxx(isnan(Sxx)) = 0;
Syy(isnan(Syy)) = 0;
Sxy(isnan(Sxy)) = 0;

% Pricipal stresses
S1 = (Sxx+Syy)/2 + sqrt((Sxx-Syy).^2/4 + Sxy.^2);
S2 = (Sxx+Syy)/2 - sqrt((Sxx-Syy).^2/4 + Sxy.^2);
% Principal orientation
pangle = atan(2*Sxy./(Sxx-Syy)) / 2;

% Save data
save(savename,'Sxx','Syy','Sxy','pangle','S1','S2','x','y');


