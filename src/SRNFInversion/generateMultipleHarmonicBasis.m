% Generates a series of harmonics basis (L=1...15) and saves them for future
% usage.
% By Hamid Laga
% This function calls genHarmBasis of Qian
%---
function generateMultipleHarmonicBasis(res, L, outDir)
%
% Res - the desired spherical resolutions
% L - array that contains the L values we would like to generate
% outDir - where to save the generated basis

if nargin < 3
   outDir =  './HarmonicBasis/';
end

if nargin<2,
    L = floor(res(1)/2);
end

for l = 2:L,  % L:L, % 24:2:L, %
    disp(['L = ' num2str(l) ' ...']);    
    [B, Bu, Bv, ~] = genHarmBasis(res, l, 0);
    save( [ outDir 'harmonic_basis_Res' num2str(res(1)) '_L' num2str(l) '.mat'], 'B', 'Bu', 'Bv');  
end

disp('Done.');