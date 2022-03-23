function f = genSurfR3(c,B)

[ndim,nPt,numB] = size(B);

% f = zeros(ndim,nPt);
% for i1=1:numB
%     f = f + c(i1).*B(:,:,i1);
% end

Bvec = reshape(B,[ndim*nPt,numB]);
f = Bvec*c(:);
f = reshape(f,[ndim,nPt]);

end