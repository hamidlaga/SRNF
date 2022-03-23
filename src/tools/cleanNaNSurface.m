function Mnew = cleanNaNSurface(M)
%
% Cleans the NaNs if any
%

[w, h, d] = size(M);

Mnew = M;

% for k=1:d,
%     Mnew(:, :, k) = medfilt2(squeeze(M(:, :, k)));
% end

for i=1:w,
    for j=1:h,
        
       for k=1:d,
           if isnan(M(i, j, k))               
               Mnew(i, j, k) = cleanImage(squeeze(M(:, :, k)), i, j);
           end
       
       end
       
    end
end