%
% The input is a surface
%
% The third outpu F is the sphere
%
function [F2, O, F] = alignSurfaceGridWithSphere(M2n, surfaceCodeDir, currDir)

[res(1), res(2), ~] = size(M2n);

[Theta, Phi] = getSphericalGrid(res(1), 0.01);
[X, Y, Z]    = spherical_to_cart_m(Theta, Phi);
F(:,:, 1) = X;
F(:,:, 2) = Y;
F(:,:, 3) = Z;

[A1,~,~] = area_surf_closed(F);
F = F / sqrt(A1); 
    
figure(33); clf;
surface(F(:,:,1),F(:,:,2),F(:,:,3));
axis equal;
cameramenu;

cd([surfaceCodeDir 'Mex\ClosedIan']);  % 'Matlab\ClosedIan']); % 
addpath([surfaceCodeDir 'Matlab\ClosedIan']);
addpath([surfaceCodeDir 'Matlab\Closed']);

[Snorm2n,A2n,multfact2n,sqrtmultfact2n] = area_surf_closed(F);
q1 = surface_to_q(Snorm2n,sqrtmultfact2n);

% the second surface
[Snorm2n,A2n,multfact2n,sqrtmultfact2n] = area_surf_closed(M2n);
q2n = surface_to_q(Snorm2n,sqrtmultfact2n);

Ot  = optimal_rot_surf_closed(q1,q2n,Theta);
cd (currDir);

O  = Ot';
F2 = rotateGrid(M2n, O);