function [Qs, Qus, Qvs] = surfaces2srnf(Ss)
    
Qs  = zeros(size(Ss));
Qus = zeros(size(Ss));
Qvs = zeros(size(Ss));

for i=1:size(Ss, 1)
    [Qs(i, :, :,:), Qus(i, :, :,:), Qvs(i, :, :,:)]  = surface2srnf(squeeze(Ss(i, :, :, :)));    
end
