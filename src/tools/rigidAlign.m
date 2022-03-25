%
% Rigid alignment in the SRNF space
%
function [M2n, O] = rigidAlign(M1, M2, surfaceCodeDir, currDir)

[res(1), res(2), ~] = size(M1);

[Theta, Phi] = getSphericalGrid(res(1), 0.01);
[X, Y, Z]    = spherical_to_cart_m(Theta, Phi);
F(:,:, 1) = X;
F(:,:, 2) = Y;
F(:,:, 3) = Z;

% figure(33); clf;
% surface(F(:,:,1),F(:,:,2),F(:,:,3));
% axis equal;
% cameramenu;

% cd([surfaceCodeDir 'Mex/ClosedIan']);  % 
cd ([surfaceCodeDir 'Matlab/ClosedIan']); % 
addpath([surfaceCodeDir 'Matlab/ClosedIan']);
addpath([surfaceCodeDir 'Matlab/Closed']);

    
[Snorm2n,A2n,multfact2n,sqrtmultfact2n] = area_surf_closed(M1);
q1 = surface_to_q(Snorm2n,sqrtmultfact2n);

% the second surface
[Snorm2n,A2n,multfact2n,sqrtmultfact2n] = area_surf_closed(M2);
q2 = surface_to_q(Snorm2n,sqrtmultfact2n);

O  = optimal_rot_surf_closed(q1,q2,Theta);


M2n = rotate3D(M2, O);

cd (currDir);