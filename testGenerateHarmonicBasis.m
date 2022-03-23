%% generating harmonic basis
%
%  You need to call generateMultipleHarmonicBasis(res, L)
%  But you need to choose the resolution and L 
% 
% In  the example below, it generates harmonic basis for surfaces of
% resolution res = 16. It also varies L from 2 to 36. You can change it if
% you want.
%
% The generated harmonic basis will be saved in the folder specified in the variable outDir 

addpath('src/');

res = 16;
L   = 2:1:36; 
outDir ='./HarmonicBasis/';

mkdir(outDir);  % makes  the output directory if it does not exist

%% Now generate the harmonic basis - the result will be saved in outDir
generateMultipleHarmonicBasis(res, L); % here it will take some time
