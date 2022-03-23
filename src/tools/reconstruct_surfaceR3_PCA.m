function [c_new,f_new] = reconstruct_surfaceR3_PCA(f_orig, B, res, Mu)

if nargin < 4,
    Mu = zeros(3, prod(res));
end

% disp('(Surface) Reconstructing surface...')

[d,nPt,nB] = size(B);
c_new = zeros(nB, 1);

for j=1:nB
    V = squeeze(B(:,:,j));
    
    c_new(j) = (f_orig(:) - Mu(:))' * V(:);   
end

f_new = genSurfR3(c_new, B) + Mu;

end