function [Sxx, Syy, Sxy, xc, yc, um, vm] = compute_stress(x,y,tx,ty,K1,K2,h,BC)
%COMPUTE_STRESS.m COMPUTE MONOLAYER STRESSES FROM TRACTIONS
% 
% [Sxx, Syy, Sxy, xc, yc, um, vm] = compute_stress(x,y,tx,ty,K1,K2,h,BC)
%
% Compute in-plane stress tensor components from a list of grid points
% (x,y) and nodal tractions (tx,ty).
% 
% INPUTS
% 
% x, y      Nodes for computing FE mesh (vector); units: m
% 
% tx, ty    Traction components at each node (vector): units: Pa. These are
%           the tractions applied by the cell monolayer to the substrate.
%           If you input tractions applied by the substrate to the cells,
%           the results will have the wrong sign.
% 
% K1, K2    Constitutive properties for monolayer. Can be either scalar 
%           (indicating constant value across the monolayer) or a vector 
%           specifying properties at each node.
%           
%              If monolayer is assumed to be a viscous fluid, these are the
%              bulk (K1) and shear (K2) viscosities.
% 
%              If monolayer is assumed to be an elastic solid, these are 
%              the bulk (K1) and shear (K2) moduli. 
% 
%              By the correspondence principle of viscoelasticity, this
%              implementation also applies for a viscoelastic fluid or
%              a viscoelastic solid.
% 
%           Units of K1 and K2 must be the same, but the units chosen will
%           not affect the results, because the results depend on a 
%           dimensionless ratio of K1 and K2, but not on K1 and K2
%           directly.
%
%           The implementation for both solid and fluid assumes no stress
%           in out-of-plane direction. However, there can be nonzero motion
%           in the out-of-plane direction, i.e., there can be nonzero
%           out-of-plane strains/strain rates. This is consistent with
%           experimental data showing cells in a monolayer change height
%           during migration.
% 
%           In both cases, the solution depends on a dimensionless ratio of 
%           K1 and K2. Therefore, the magnitudes of K1 and K2 don't affect
%           the results; only the dimensionless ratio does.
%           
% h         Height of monolayer. Can be either a scalar or a vector
%           specifing the height at each node. Units: m
% 
% BC        Boundary conditions to apply. There are several options:
%               'left-right' - Leftmost and rightmost points are 
%                              constrained, giving 3 BCs to account for 3 
%                              degrees of freedom
%               'top' - top boundary is onstrained
%               'bottom' - bottom boundary is onstrained
%               'left' - left boundary is onstrained
%               'right' - right boundary is onstrained
%               'all' - The boundary of the monolayer is onstrained (not
%                       recommended)
%               'none' - No boundary conditions (useful for testing, but
%                        not recommended)
%          Note that the results are very sensitive the the boundary
%          conditions. Typically, when the code is run on a cell island,
%          the BC 'left-right' gives a stress tensor having a positive
%          trace (implying tension) at all locations in the monolayer; by
%          contrast, the BC 'all' gives a stress tensor whose trace
%          averages to 0 across the monolayer. Therefore the option 'all'
%          is not recommended. The options 'top', 'bottom', 'left', and
%          'right' are intended to be used for cell monolayers that are not
%          islands but rather occupy the top, bottom, left, or right side
%          of an image. In this case, the top, bottom, left, or right is 
%          referred to as an "optical boundary" by Tambe et al, Plos One,
%          2013, e55172. See this reference for details on sensitivity of
%          the results to the presence of an optical boundary.
%
% OUTPUTS
% 
% Sxx, Syy, Sxy     Stress tensor components at element centers (Pa)
% xc, yc            Positions of element centers (m)
% um, vm            Vector field at nodes
%                   If K1 & K2 have units of Pa-s, then um, vm are m/s
%                   If K1 & K2 have units of Pa, then units of um, vm are m
%                   
% 
% NOTES
% 
% Inputs are all vectors. This script computes stresses using every input 
% point, input only the points that correspond to the cell monolayer. The 
% code assumes x and y have constant and equal spacing.
% 
% Outputs are vectors of the same size as the inputs. Stresses are computed
% at the center of each element, (xc, yc). The components of the vector 
% field (um, vm) are computed at the nodes.
% 
% This function uses a Delaunay triangulation of the input coordinates
% (x,y) to create 3-node triangular elements.
% 
% This code follows the procedure of the following papers:
%  Tambe et al, Nat Mater, 2011, 10:469-475
%  Tambe et al, Plos One, 2013, e55172
% 
% Written by Jacob Notbohm and Aashrith Saraswathibhatla, University of 
% Wisconsin-Madison, 2015-2019
% 

% This script is written generally to account for different height and
% constitutive properties across the cell monolayer.

%% --- PRELIMINARIES ---

tic

% Compute grid spacing. This requires a uniform grid everywhere with equal
% spacing in x and y directions
xx = unique(sort(x));
grid_spacing = xx(2)-xx(1);

% Turn possible scalar inputs into vectors
if numel(K1)==1
    K1 = K1*ones(size(x));
end
if numel(K2)==1
    K2 = K2*ones(size(x));
end
if numel(h)==1 % Height
    h = h*ones(size(x));
end

% Delaunay triangulation. This outputs a Nx3 array TRI where N is the
% number of elements. The contents of TRI are the indices the vector inputs
% x and y that give the coordinates of the 3 nodes of the n-th element.
DT = delaunayTriangulation(x,y);
TRI = DT.ConnectivityList;
% Note: It's important to avoid long, skinny triangles. There's a check for 
% this a few lines below. 

%% --- PROPERTIES OF ELEMENTS ---

% x,y node positions of elements
x1 = x(TRI(:,1));
y1 = y(TRI(:,1));
x2 = x(TRI(:,2));
y2 = y(TRI(:,2));
x3 = x(TRI(:,3));
y3 = y(TRI(:,3));

% Element areas
A_el = ( (x2.*y3-x3.*y2) + (x3.*y1-x1.*y3) + (x1.*y2-x2.*y1) )/2 ;

% Check for small areas resulting from bad triangulation. Elements with 
% very small area are normally on the edges of the domain. Thus, they can 
% be removed without any problems.
idx = A_el < 1e-6*median(A_el) ;
if any(idx)
    % Remove from node list and x,y positions
    TRI = TRI(~idx,:);
    x1 = x(TRI(:,1));
    y1 = y(TRI(:,1));
    x2 = x(TRI(:,2));
    y2 = y(TRI(:,2));
    x3 = x(TRI(:,3));
    y3 = y(TRI(:,3));
    A_el = ( (x2.*y3-x3.*y2) + (x3.*y1-x1.*y3) + (x1.*y2-x2.*y1) )/2 ;
end

% Define element constitutive props to be mean of those at nodes.
K1_el = ( K1(TRI(:,1)) + K1(TRI(:,2)) + K1(TRI(:,3)) )/3 ;
K2_el = ( K2(TRI(:,1)) + K2(TRI(:,2)) + K2(TRI(:,3)) )/3 ;
% Element heights. Use mean value of heights at nodes. This is approximate
% if constitutive properties and h vary over space
h_el = ( h(TRI(:,1)) + h(TRI(:,2)) + h(TRI(:,3)) )/3 ;


%% --- ASSEMBLE STIFFNESS MATRIX ---

% Preallocate vectors for stiffness matrix. Vectors will be used to
% assemble a sparse array.
Nel = length(TRI(:,1)); % Number of elements
idx_row = zeros(1,Nel*36); % rows of nonzero elements of stiffness matrix
idx_col = zeros(1,Nel*36); % cols of nonzero elements of stiffness matrix
Kv = zeros(1,Nel*36); % nonzero elements of stiffness matrix

for k=1:Nel
    % Indices of nodes
    n1 = TRI(k,1);
    n2 = TRI(k,2);
    n3 = TRI(k,3);
    % Assemble matrix of consititutive properties. This follows the 
    % implementation of plane stress, in which stresses in the out-of-plane
    % direction are zero.
    E_el_matrix = 2*K2_el(k)/(3*K1_el(k) + 4*K2_el(k))*...
                    [ 6*K1_el(k)+2*K2_el(k), 3*K1_el(k)-2*K2_el(k), 0 ; ...
                      3*K1_el(k)-2*K2_el(k), 6*K1_el(k)+2*K2_el(k), 0 ; ...
                      0,                     0,                     3*K1_el(k)+4*K2_el(k) ];
    % Matrix giving the gradient of the vector field: epsilon_el = B*u_el
    B = 1/(2*A_el(k)) * [ y2(k)-y3(k), 0, y3(k)-y1(k), 0, y1(k)-y2(k), 0 ; ...
                    0, x3(k)-x2(k), 0, x1(k)-x3(k), 0, x2(k)-x1(k) ; ...
                    x3(k)-x2(k), y2(k)-y3(k), x1(k)-x3(k), y3(k)-y1(k), x2(k)-x1(k), y1(k)-y2(k) ];
    % Stiffness matrix for k-th element
    K_k = A_el(k)*h_el(k)*B'*E_el_matrix*B;
    
    % % Test element stiffness matrix by computing e-values. 3 e-values must
    % % be positive; 3 must be zero
    % lambda = eig(K_k)
    
    % This stiffness matrix is 6x6 (3 nodes, 2 dof each)
    % Corresponding vector field is [ux1, uy1, ux2, uy2, ux3, uy3]'
    % To add to global stiffness matrix, forces correspond to rows;
    % vector field corresponds to columns
    
    % ROW 1
    idx_row(36*(k-1)+1) = 2*n1-1; 
    idx_col(36*(k-1)+1) = 2*n1-1;
    Kv(36*(k-1)+1) = K_k(1,1);
    idx_row(36*(k-1)+2) = 2*n1-1; 
    idx_col(36*(k-1)+2) = 2*n1;
    Kv(36*(k-1)+2) = K_k(1,2);
    idx_row(36*(k-1)+3) = 2*n1-1;
    idx_col(36*(k-1)+3) = 2*n2-1;
    Kv(36*(k-1)+3) = K_k(1,3);
    idx_row(36*(k-1)+4) = 2*n1-1;
    idx_col(36*(k-1)+4) = 2*n2;
    Kv(36*(k-1)+4) = K_k(1,4);
    idx_row(36*(k-1)+5) = 2*n1-1;
    idx_col(36*(k-1)+5) = 2*n3-1;
    Kv(36*(k-1)+5) = K_k(1,5);
    idx_row(36*(k-1)+6) = 2*n1-1;
    idx_col(36*(k-1)+6) = 2*n3;
    Kv(36*(k-1)+6) = K_k(1,6);
    
    % ROW 2
    idx_row(36*(k-1)+7) = 2*n1; 
    idx_col(36*(k-1)+7) = 2*n1-1;
    Kv(36*(k-1)+7) = K_k(2,1);
    idx_row(36*(k-1)+8) = 2*n1; 
    idx_col(36*(k-1)+8) = 2*n1;
    Kv(36*(k-1)+8) = K_k(2,2);
    idx_row(36*(k-1)+9) = 2*n1;
    idx_col(36*(k-1)+9) = 2*n2-1;
    Kv(36*(k-1)+9) = K_k(2,3);
    idx_row(36*(k-1)+10) = 2*n1;
    idx_col(36*(k-1)+10) = 2*n2;
    Kv(36*(k-1)+10) = K_k(2,4);
    idx_row(36*(k-1)+11) = 2*n1;
    idx_col(36*(k-1)+11) = 2*n3-1;
    Kv(36*(k-1)+11) = K_k(2,5);
    idx_row(36*(k-1)+12) = 2*n1;
    idx_col(36*(k-1)+12) = 2*n3;
    Kv(36*(k-1)+12) = K_k(2,6);
    
    % ROW 3
    idx_row(36*(k-1)+13) = 2*n2-1;
    idx_col(36*(k-1)+13) = 2*n1-1;
    Kv(36*(k-1)+13) = K_k(3,1);
    idx_row(36*(k-1)+14) = 2*n2-1; 
    idx_col(36*(k-1)+14) = 2*n1;
    Kv(36*(k-1)+14) = K_k(3,2);
    idx_row(36*(k-1)+15) = 2*n2-1;
    idx_col(36*(k-1)+15) = 2*n2-1;
    Kv(36*(k-1)+15) = K_k(3,3);
    idx_row(36*(k-1)+16) = 2*n2-1;
    idx_col(36*(k-1)+16) = 2*n2;
    Kv(36*(k-1)+16) = K_k(3,4);
    idx_row(36*(k-1)+17) = 2*n2-1;
    idx_col(36*(k-1)+17) = 2*n3-1;
    Kv(36*(k-1)+17) = K_k(3,5);
    idx_row(36*(k-1)+18) = 2*n2-1;
    idx_col(36*(k-1)+18) = 2*n3;
    Kv(36*(k-1)+18) = K_k(3,6);
    
    % ROW 4
    idx_row(36*(k-1)+19) = 2*n2;
    idx_col(36*(k-1)+19) = 2*n1-1;
    Kv(36*(k-1)+19) = K_k(4,1);
    idx_row(36*(k-1)+20) = 2*n2; 
    idx_col(36*(k-1)+20) = 2*n1;
    Kv(36*(k-1)+20) = K_k(4,2);
    idx_row(36*(k-1)+21) = 2*n2;
    idx_col(36*(k-1)+21) = 2*n2-1;
    Kv(36*(k-1)+21) = K_k(4,3);
    idx_row(36*(k-1)+22) = 2*n2;
    idx_col(36*(k-1)+22) = 2*n2;
    Kv(36*(k-1)+22) = K_k(4,4);
    idx_row(36*(k-1)+23) = 2*n2;
    idx_col(36*(k-1)+23) = 2*n3-1;
    Kv(36*(k-1)+23) = K_k(4,5);
    idx_row(36*(k-1)+24) = 2*n2;
    idx_col(36*(k-1)+24) = 2*n3;
    Kv(36*(k-1)+24) = K_k(4,6);
    
    % ROW 5
    idx_row(36*(k-1)+25) = 2*n3-1;
    idx_col(36*(k-1)+25) = 2*n1-1;
    Kv(36*(k-1)+25) = K_k(5,1);
    idx_row(36*(k-1)+26) = 2*n3-1; 
    idx_col(36*(k-1)+26) = 2*n1;
    Kv(36*(k-1)+26) = K_k(5,2);
    idx_row(36*(k-1)+27) = 2*n3-1;
    idx_col(36*(k-1)+27) = 2*n2-1;
    Kv(36*(k-1)+27) = K_k(5,3);
    idx_row(36*(k-1)+28) = 2*n3-1;
    idx_col(36*(k-1)+28) = 2*n2;
    Kv(36*(k-1)+28) = K_k(5,4);
    idx_row(36*(k-1)+29) = 2*n3-1;
    idx_col(36*(k-1)+29) = 2*n3-1;
    Kv(36*(k-1)+29) = K_k(5,5);
    idx_row(36*(k-1)+30) = 2*n3-1;
    idx_col(36*(k-1)+30) = 2*n3;
    Kv(36*(k-1)+30) = K_k(5,6);
    
    % ROW 6
    idx_row(36*(k-1)+31) = 2*n3;
    idx_col(36*(k-1)+31) = 2*n1-1;
    Kv(36*(k-1)+31) = K_k(6,1);
    idx_row(36*(k-1)+32) = 2*n3; 
    idx_col(36*(k-1)+32) = 2*n1;
    Kv(36*(k-1)+32) = K_k(6,2);
    idx_row(36*(k-1)+33) = 2*n3;
    idx_col(36*(k-1)+33) = 2*n2-1;
    Kv(36*(k-1)+33) = K_k(6,3);
    idx_row(36*(k-1)+34) = 2*n3;
    idx_col(36*(k-1)+34) = 2*n2;
    Kv(36*(k-1)+34) = K_k(6,4);
    idx_row(36*(k-1)+35) = 2*n3;
    idx_col(36*(k-1)+35) = 2*n3-1;
    Kv(36*(k-1)+35) = K_k(6,5);
    idx_row(36*(k-1)+36) = 2*n3;
    idx_col(36*(k-1)+36) = 2*n3;
    Kv(36*(k-1)+36) = K_k(6,6);
    
end

% Create stiffness matrix as a sparse array
K = sparse(idx_row, idx_col, Kv);

%% --- ASSEMBLE FORCE VECTOR ---

N = length(x); % N is number of nodes

% Note: Do not scale tractions by height here, because height is accounted
% for in the stiffness matrix.
b = zeros(2*N,1);
b(1:2:2*N-1) = tx;
b(2:2:2*N) = ty;

% Scale tractions by area to convert from stress to force.
area = (grid_spacing)^2;
b = b*area;

% The tractions that have been input into this script are the tractions
% applied by the cells to the substrate. We want the tractions applied by
% the substrate to the cells.
b = -b;

%% --- APPLY BOUNDARY CONDITIONS ---

% NOTE: Solution is not unique:
% 1. det(K)=0
% 2. Null space is made of 3 unique vectors. Thus, we need to specify 3
% dofs (2 for translation, 1 for rotation).
% 3. 3 e-values are 0. Rest of e-values are postitive. This is correct.
% Thus, 3 boundary conditions are required.

% Various boundary conditions are implemented here, allowing the user to
% consider their effects in the region of interest of the cell monolayer

switch BC
    case 'none'
    % --- NO BOUNDARY CONDITION (START) ---
    % For the island geometry, there can be rotations which cause large shears
    % to appear. It's good to check that the rotations are small by solving
    % with no boundary conditions. It's not recommended to use this case in
    % general, as the solution is not unique.
    
    % This only requires updating the notation of the stiffness and force
    % arrays
    K_BC = K;
    b_BC = b;
    % --- NO BOUNDARY CONDITION (END) ---
    
    case 'left-right'
    % --- BOUNDARY CONDITION 'LEFT-RIGHT' (START) ---
    % For the island geometry, it's common to use pins on the left and
    % right. Do this by adding 3 Lagrange multipliers
    K_BC = sparse(2*N+3,2*N+3);
    K_BC(1:2*N,1:2*N) = K_BC(1:2*N,1:2*N) + K;
    % Get node numbers of left and right sides
    [~, n_left] = min(x);
    [~, n_right] = max(x);
    
    % If there exist more than 1 node on the left or right, choose the one
    % corresponding to the mean y position.
    % Begin with left side
    idx = find(abs(x-x(n_left))/grid_spacing < 1e-6);
    if length(idx)>1
        ymean = mean(y(idx));
        [~, i2] = min(abs(y(idx)-ymean));
        n_left = idx(i2);
    end
    % Repeat for right side
    idx = find(abs(x-x(n_right))/grid_spacing < 1e-6);
    if length(idx)>1
        ymean = mean(y(idx));
        [~, i2] = min(abs(y(idx)-ymean));
        n_right = idx(i2);
    end
    
    % Constrain dof 2 on left node
    K_BC(2*N+1,2*n_left) = 1;
    K_BC(2*n_left,2*N+1) = 1;
    % Constrain dof 2 on right node
    K_BC(2*N+2,2*n_right) = 1;
    K_BC(2*n_right,2*N+2) = 1;
    % Constrain dof 1 on right node
    K_BC(2*N+3,2*n_right-1) = 1;
    K_BC(2*n_right-1,2*N+3) = 1;
    
    % Create force array including zero forces on the 3 lagrange multipliers
    b_BC = zeros(2*N+3,1);
    b_BC(1:2*N) = b;
    % --- BOUNDARY CONDITION 'LEFT-RIGHT' (END)----
    
    case 'all'
    % --- BOUNDARY CONDITION 'ALL' (START) ---
    % In this boundary condition, we pin all the boundary points of the
    % cell monolayer
    
    % Identify nodes on boundary. The code below uses Matlab scripts to
    % find the boundary of the triangulation. Note that this only
    % identifies nodes on the convex hull
    boundarynodes = freeBoundary(DT);
    boundarynodes = unique(boundarynodes(:));
    
    K_BC = sparse(2*N+2*length(boundarynodes), 2*N+2*length(boundarynodes));
    K_BC(1:2*N,1:2*N) = K_BC(1:2*N,1:2*N) + K;
    for k=1:length(boundarynodes)
        K_BC(2*N+2*k-1, 2*boundarynodes(k)-1) = 1;  % dof 1
        K_BC(2*boundarynodes(k)-1, 2*N+2*k-1) = 1;
        K_BC(2*N+k, 2*boundarynodes(k)) = 1;        % dof 2
        K_BC(2*boundarynodes(k), 2*N+k) = 1;
    end
    
    % Create force array including zero forces on the Lagrange multipliers
    b_BC = zeros(2*N+2*length(boundarynodes), 1);
    b_BC(1:2*N) = b; 
    
    % An alternative method is below. This identifies more boundary nodes
    % in some cases.
%     xu = unique(x);
%     [y_left,~] = find(abs(x-xu(1))/grid_spacing < 1); % finding the y points of the left x
%     [y_right,~] = find(abs(x-xu(end))/grid_spacing < 1); % finding the y points of the right x
%     [~,n_right] = max(x);
%     K_BC  = sparse(2*N+4*(length(xu)-2)+2*length(y_left)+2*length(y_right),...
%         2*N+4*(length(xu)-2)+2*length(y_left)+2*length(y_right));
%     K_BC(1:2*N,1:2*N) = K_BC(1:2*N,1:2*N) + K;
%     for i = 2:length(xu)-1
%         row = find(abs(x-xu(i))/grid_spacing < 1);
%         ymin = min(row);
%         ymax = max(row);
%         % dof of the ymin at a certain x
%         K_BC(2*N+i-1,2*ymin) = 1; % dof 2
%         K_BC(2*ymin,2*N+i-1) = 1; % dof 2
%         K_BC(2*N+length(xu)+i-3,2*ymin-1)=1;
%         K_BC(2*ymin-1,2*N+length(xu)+i-3)=1;
%         % dof of the ymax at a certain x
%         K_BC(2*N+2*length(xu)+i-5,2*ymax) = 1;
%         K_BC(2*ymax,2*N+2*length(xu)+i-5) = 1;
%         K_BC(2*N+3*length(xu)+i-7,2*ymax-1) = 1;
%         K_BC(2*ymax-1,2*N+3*length(xu)+i-7) = 1;
%         clear row
%     end
%     for i=1:length(y_left)
%         % pinning the left edge
%         K_BC(2*N+4*(length(xu)-2)+i,2*i) = 1;
%         K_BC(2*i,2*N+4*(length(xu)-2)+i) = 1;
%         K_BC(2*N+4*(length(xu)-2)+length(y_left)+i,2*i-1) = 1;
%         K_BC(2*i-1,2*N+4*(length(xu)-2)+length(y_left)+i) = 1;
%     end
%     for i=1:length(y_right)
%         % pinning the right edge
%         K_BC(2*N+4*(length(xu)-2)+2*length(y_left)+i,2*i+2*n_right) = 1;
%         K_BC(2*i+2*n_right,2*N+4*(length(xu)-2)+i+2*length(y_left)) = 1;
%         K_BC(2*N+4*(length(xu)-2)+length(y_right)+2*length(y_left)+i,2*i-1+2*n_right) = 1;
%         K_BC(2*i-1+2*n_right,2*N+4*(length(xu)-2)+length(y_right)+2*length(y_left)+i) = 1;
%     end
%     
%     % Create force array including zero forces on the many lagrange multipliers
%     b_BC = zeros(2*N+4*(length(xu)-2)+2*length(y_left)+2*length(y_right),1);
%     b_BC(1:2*N) = b;
    % --- BOUNDARY CONDITION 'ALL' (END) ---
    
    case 'left'
        % --- BOUNDARY CONDITION 'LEFT' (START) ---
        % pinning the left boundary
        xu = unique(x);
        row = find(abs(x-xu(1))/grid_spacing < 1);
        K_BC = sparse(2*N+2*length(row),2*N+2*length(row));
        K_BC(1:2*N,1:2*N) = K_BC(1:2*N,1:2*N) + K;
        for i=1:length(row)
            K_BC(2*N+i,2*i) = 1;
            K_BC(2*i,2*N+i) = 1;
            K_BC(2*N+length(row)+i,2*i-1) = 1;
            K_BC(2*i-1,2*N+length(row)+i) = 1;
        end
        
        % Force array with lagrange multipliers
        b_BC = zeros(2*N+2*length(row),1);
        b_BC(1:2*N) = b;
        % --- BOUNDARY CONDITION 'LEFT' (END) ---
        
    case 'right'
        % --- BOUNDARY CONDITION 'RIGHT' (START) ---
        % pinning the right boudary
        xu = unique(x);
        row = find(abs(x-xu(end))/grid_spacing < 1);
        [~,n_right] = max(x);
        K_BC = sparse(2*N+2*length(row),2*N+2*length(row));
        K_BC(1:2*N,1:2*N) = K_BC(1:2*N,1:2*N) + K;
        for i=1:length(row)
            K_BC(2*N+i,2*n_right+2*i) = 1;
            K_BC(2*n_right+2*i,2*N+i) = 1;
            K_BC(2*N+i+length(row),2*n_right+2*i-1) = 1;
            K_BC(2*n_right+2*i-1,2*N+i+length(row)) = 1;
        end 
        
        % Force array with lagrange multipliers
        b_BC = zeros(2*N+2*length(row),1);
        b_BC(1:2*N) = b;
        % --- BOUNDARY CONDITION 'RIGHT' (END) ---
        
    case 'top'
        % --- BOUNDARY CONDITION 'TOP' (START) ---
        % pinning the top boundary
        ymax = max(y);
        row = find(abs(y-ymax)/grid_spacing < 1);
        K_BC = sparse(2*N+2*length(row),2*N+2*length(row));
        K_BC(1:2*N,1:2*N) = K_BC(1:2*N,1:2*N) + K;
        for i=1:length(row)
            K_BC(2*N+i,2*row(i)) = 1;
            K_BC(2*row(i),2*N+i) = 1;
            K_BC(2*N+length(row)+i,2*row(i)-1) = 1;
            K_BC(2*row(i)-1,2*N+length(row)+i) = 1;
        end
        
        % Force array with lagrange multipliers
        b_BC = zeros(2*N+2*length(row),1);
        b_BC(1:2*N) = b;
        % --- BOUNDARY CONDITION 'TOP' (END) ---
        
    case 'bottom'
        % --- BOUNDARY CONDITION 'BOTTOM' (START) ---
        % pinning the bottom boundary
        ymin = min(y);
        row = find(abs(y-ymin)/grid_spacing < 1);
        K_BC = sparse(2*N+2*length(row),2*N+2*length(row));
        K_BC(1:2*N,1:2*N) = K_BC(1:2*N,1:2*N) + K;
        for i=1:length(row)
            K_BC(2*N+i,2*row(i)) = 1;
            K_BC(2*row(i),2*N+i) = 1;
            K_BC(2*N+length(row)+i,2*row(i)-1) = 1;
            K_BC(2*row(i)-1,2*N+length(row)+i) = 1;
        end
        
        % Force array with lagrange multipliers
        b_BC = zeros(2*N+2*length(row),1);
        b_BC(1:2*N) = b;
        % --- BOUNDARY CONDITION 'BOTTOM' (END) ---
end


%% --- SOLVE MATRIX EQUATION ---

% Solve
u = K_BC\b_BC;

% Keep only the values that aren't associated with a lagrange multiplier
um = u(1:2:2*N-1);
vm = u(2:2:2*N);


%% --- COMPUTE STRESSES IN EACH ELEMENT ---
% Note: The code below computes the element stresses with a simple average,
% which isn't the most elegant method. A better stress recovery procedure 
% may be to compute nodal point stresses using an averaging procedure.

% Preallocate arrays for stress components
Sxx = zeros(Nel,1);
Syy = Sxx;
Sxy = Sxx;

for k=1:Nel
    % Indices of nodes
    n1 = TRI(k,1);
    n2 = TRI(k,2);
    n3 = TRI(k,3);
    % Vector field of 3 nodes
    u_k = [ u(2*n1-1) ; ...
            u(2*n1); ...
            u(2*n2-1) ; ...
            u(2*n2); ...
            u(2*n3-1) ; ...
            u(2*n3) ] ;
    % Assemble matrix of constitutive properties
    E_el_matrix = 2*K2_el(k)/(3*K1_el(k) + 4*K2_el(k))*...
                    [ 6*K1_el(k)+2*K2_el(k), 3*K1_el(k)-2*K2_el(k), 0 ; ...
                      3*K1_el(k)-2*K2_el(k), 6*K1_el(k)+2*K2_el(k), 0 ; ...
                      0,                     0,                     3*K1_el(k)+4*K2_el(k) ];
    % Matrix relating vector field to its gradient
    B = 1/(2*A_el(k)) * [ y2(k)-y3(k), 0, y3(k)-y1(k), 0, y1(k)-y2(k), 0 ; ...
                    0, x3(k)-x2(k), 0, x1(k)-x3(k), 0, x2(k)-x1(k) ; ...
                    x3(k)-x2(k), y2(k)-y3(k), x1(k)-x3(k), y3(k)-y1(k), x2(k)-x1(k), y1(k)-y2(k) ];
	
	% Stress components
    S_k = E_el_matrix*B*u_k;
    
    % Add to stress arrays
    Sxx(k) = S_k(1);
    Syy(k) = S_k(2);
    Sxy(k) = S_k(3);
    
end

% Compute centers of elements
xc = (x1+x2+x3)/3;
yc = (y1+y2+y3)/3;


tm = toc;






