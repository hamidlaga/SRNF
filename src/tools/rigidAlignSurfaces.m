function Snew = rigidAlignSurfaces(Source, Target)
% Rigid aligns source Source to target Target
% Assumes that the two surfaces are of the size size nxnx3

[n,m,~] = size(Target);
P1 = reshape(Target, n*m, []);
P1 = P1';

Q1 = reshape(Source, n*m, []); 
Q1 = Q1';

A  = P1 * Q1';
[U, S, V] = svd(A); 
O         = U*V';         % aligns Q on P

Snew      = rotate3D(Source, O);