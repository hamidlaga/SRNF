function Xmean = computeSurfaceMean(MM)

Xmean = [];

if isempty(MM)
   return;
end

[nsurfaces, h, w, dim] = size(MM);
if dim<3
    return;
end

Xmean = squeeze(MM(1, :, :, :));
if nsurfaces <2,
    return;
end

QQ = zeros(nsurfaces, h, w, 3);
QQu = zeros(nsurfaces, h, w, 3);
QQv = zeros(nsurfaces, h, w, 3);

for i=1:nsurfaces
   [QQ(i, : ,:, :),  QQu(i, : ,:, :),  QQv(i, : ,:, :)] = surface2srnf(squeeze(MM(i, :, :, :) ));
end

%% Computing the mean using an iterative approach
[meanQ, meanQu, meanQv] = computSRNFmean(QQ, QQu, QQv);

%% Inversion
for i=1:size(QQ, 1)
    Xmean = invertSRNF_by_SRVF_inversion(meanQ, meanQu, meanQv);
end

%% Visualization
% figure(100);
% hold on;
% for i=1:size(SS, 1)
%     S = squeeze(SS(i, :, :, :));
%     surface(squeeze(S(:,:, 1)) + (i-1)*1.5, squeeze(S(:,:, 2)), squeeze(S(:,:, 3))); %  , 'CDATA', K);   
% end
% campos([0.6450    7.4510   13.3250]);
% axis equal;axis off;
% cameramenu;
