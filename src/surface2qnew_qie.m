%
% Converts the surface into qmap . uses Xie's representation
%
function q  = surface2qnew_qie(f, res, mysmall)

if nargin < 3,
    mysmall = eps;
end

[Theta,~] = genGridSphr(res);

sinTheta  = genSinPhi(Theta, mysmall);
N       = prod(res);
d       = 3;

[Mu, Mv] = partialF(f, sinTheta, d, N, res);

[q,  ~, ~, ~] = map2qnew(Mu, Mv, N);
