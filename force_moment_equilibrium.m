function [tx_equil, ty_equil] = force_moment_equilibrium(tx,ty,r,theta)
% 
% [tx_equil, ty_equil] = force_moment_equilibrium(tx,ty,r,theta)
% 
% Adjust x and y tractions to approach force and moment equilibrium.
% Run this function iteratively for best results.
% 
% INPUTS
%   tx, ty: x and y components of tractions
%   r, theta: radial and angular position of traction components
% 
% OUTPUTS
%   tx_equil, ty_equil: force and moment equilibrated x and y tractions
% 
% 
% Written by Jacob Notbohm, University of Wisconsin-Madison 2015-2020
% 


% Convert to polar coordinates
tr = tx.*cos(theta) + ty.*sin(theta);
ttheta = -tx.*sin(theta) + ty.*cos(theta);

% Net moment
M = sum( ttheta(~isnan(ttheta)).*r(~isnan(ttheta)) );

% Number of non-nan points used in moment sum
N = length(find(~isnan(ttheta)));

% Correct moments with a r-weighted correction
ttheta_mzero = ttheta - (M/N)./r;

% Convert back to x and y coordinates
tx_mzero = tr.*cos(theta) - ttheta_mzero.*sin(theta);
ty_mzero = tr.*sin(theta) + ttheta_mzero.*cos(theta);

% Correct for force balance by subtracting off mean
tx_equil = tx_mzero - mean(tx_mzero(~isnan(tx_mzero)));
ty_equil = ty_mzero - mean(ty_mzero(~isnan(ty_mzero)));
