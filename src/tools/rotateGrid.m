%
% Rotate the grid of F1 with a rotation O
%
%
function F1n = rotateGrid(F1, O)

% Grid for computing the deformation
[N, M, ~] = size(F1);
res       = [N, M];

[thmat, phimat] = genGridSphr(res, 0.01);

[Xp, Yp, Zp] = spherical_to_cart_m(thmat, phimat);

Grid(:,:,1) = Xp;
Grid(:,:,2) = Yp;
Grid(:,:,3) = Zp;

Gridn = rotate3D(Grid, O');

Xp = squeeze(Gridn(:,:,1));
Yp = squeeze(Gridn(:,:,2));
Zp = squeeze(Gridn(:,:,3));

[Th,Ph,~] = cartesian_to_sph_m(Xp, Yp, Zp);
gamnew(:, :, 1) = Th;
gamnew(:, :, 2) = Ph;
     
F1n = Apply_Gamma_Surf_Closed_new(F1, thmat, phimat, gamnew);


%%
function Fnew = Apply_Gamma_Surf_Closed_new(F,Theta,Phi,gam)

Fnew(:,:,1) = interp2(Theta,Phi,F(:,:,1),gam(:,:,1),gam(:,:,2),'spline');
Fnew(:,:,2) = interp2(Theta,Phi,F(:,:,2),gam(:,:,1),gam(:,:,2),'spline');
Fnew(:,:,3) = interp2(Theta,Phi,F(:,:,3),gam(:,:,1),gam(:,:,2),'spline');

% %%
% function qnew=rotate3Dn(q,Ot)
% 
% [n,t,d]=size(q);
% 
% for i=1:n
%     for j=1:t
%         qnew(i,j,:)=Ot*squeeze(q(i,j,:));
%     end
% end