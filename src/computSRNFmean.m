function [meanQ, meanQu, meanQv] = computSRNFmean(QQ, QQu, QQv)

nSurfaces = size(QQ, 1);

meanQ = squeeze(QQ(1, :, :, :));
meanQu = squeeze(QQu(1, :, :, :));
meanQv = squeeze(QQv(1, :, :, :));

for i=2:nSurfaces
    
    sc = [0, 1.0/i, 1];    
    
    [meanq, meanqu, meanqv] = computeGeodesicSRNF(meanQ, meanQu, meanQv, ....
                                                  squeeze(QQ(i, :, :,:)), squeeze(QQu(i, :, :,:)), squeeze(QQv(i, :, :,:)), ...
                                                  sc);
    meanQ  = squeeze(meanq(2, :, :, :));
    meanQu = squeeze(meanqu(2, :, :, :));
    meanQv = squeeze(meanqv(2, :, :, :));                                                  
    
end