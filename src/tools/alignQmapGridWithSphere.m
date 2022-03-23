%
% The input is a the qmap of a surface
%
function [q2n, O] = alignQmapGridWithSphere(q2, q1)

eps = 0.01;

[res(1), res(2), ~] = size(q2);

% Canonical spherical grid
[Theta, Phi] = getSphericalGrid(res(1), eps);
[X, Y, Z]    = spherical_to_cart_m(Theta, Phi);
F(:,:, 1) = X;
F(:,:, 2) = Y;
F(:,:, 3) = Z;

% Qmap of a sphere    
fn        = mySurf(F, 0);
sinTheta  = genSinPhi(Theta, eps);
N         = prod(res);
d         = 3;
    
if nargin < 2,
    [Mu, Mv] = partialF(fn, sinTheta, d, N, res);
    [qf,  ~, ~, ~] = map2qnew(Mu, Mv, N);

    q1  = reshape(qf, [3, res(1), res(2)]);
    q1  = permute(q1,[2,3,1]);    
end

% alignment
% size(q1)
% size(q2)

Ot  = optimal_rot_surf_closed(q1,q2,Theta);

O   = Ot';
q2n = rotateGrid(q2, O);