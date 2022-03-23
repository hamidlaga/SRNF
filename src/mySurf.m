function [f,res] = mySurf(f,flip_param)
%INPUTS:
% f: array of a single surface, dimension could be [d,n1,n2] or [n1,n2,d]
% flip_param: bool, if TRUE, then the parametrization is flipped, need to
% fix
%OUPUTS:
% f: new array of surface to analyze, dim(f) = [d,ncol*nrow]
% res: stores dim(f) = [nrow,ncol]
%USAGE: 
% [f,res] = mySurf(f);
% [f,res] = mySurf(f,1);

d = 3;
nF = size(f);

if (nF(3) == d)
%     f = closeSurface(f);
elseif (nF(1) == d)
    f = permute(f,[2,3,1]);
%     f = closeSurface(f);
else
    error('Wrong dimension of surface.')
end

f = permute(f,[3,1,2]);
nF = size(f);

if (nargin == 2) && (flip_param)
    f = permute(f,[1,3,2]);
end

show_fig = 0;
if (show_fig)
    figure; dispSurf(f(1,:,:),f(2,:,:),f(3,:,:));
    cameramenu
end

res = [nF(2),nF(3)];
f = reshape(f,[d,prod(res)]);

% disp('Surface dimension standardized: dim(f) = [d,n1*n2]');

end