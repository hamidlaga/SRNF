%% README
% This file demonstrate how to use the SRNF inversion procedure.
%
% In this example, we show how to compute the invert of an SRNF map when an
% an initial surface is not available but the space of such shapes can be represented with PCA basis.
% This is quite common for specific classes of shapes such as human bodies and faces where a morphable model can be learned from exemplars.
%
% In this example, we focus on human body shapes for  PCA basis are
% available
%
% 
% REQUIREMENTS
% - PCA Basis for the class of objects you would like to analyze. For instance,
%   if dealing with body shapes, you need to create PCA basis of humans 
%   (by performing PCA on a dataset of human body shapes)
% 
%   The script for this is testGeneratePCABasis.m. You need to edit it and
%   update the paths to point to your dataset.
%
%   In this distribution, we include PCA basis learned from the DFAUSt human body shapes:
%       https://dfaust.is.tue.mpg.de/
%
%   If you are using a different class of shapes, or a different cohort of
%   body shapes, you will need to generate the corresponding PCA basis
% 
% - Compile the mex files inside the folder src/ (precompiled files are
% included for Windows (32, 64) and Mac (32, 64))
%       mex innerS2.cpp myVector.cpp
%       mex  map2qnew.cpp myVector.cpp   
%       mex  gradInvQnew.cpp myVector.cpp
%
% Preparation of the surfaces:
% - Spherical parameterization - All surfaces need to be spherically
%   parameterized. The examples provided here are already parameterized. 
%   If you want to use your owen surfaces, you need to parameterize them using the
%   code available at: https://github.com/hamidlaga/SphericalParameterization  
%
% - Multiscale representation using spherical harmonics - Not needed when
% using PCA basis for the SRNF inversion.
%
%% Example
% The script below takes the surface of a human body shape, computes its SRNF and then tries to
% recover the initial surface by SRNF inversion. It uses the mean surface
% (provided in the PCA model) as initialisation
%
% Note that in the example below you may need to edit the paths depending
% on where you saved the code and the data.
%
% Requirements
% - PCA basis should be pre-generated for the class of shapes you are dealing with. 
%   In this distribution, we provide PCA basis learned from a subset of the
%   DFAUST data set
%
%   If you want to use a different folder, please change the variable
%     params.PCA_BASIS_FNAME_PREFIX = [params.homeDir  'PCABasis/DFAUST_PCA_with_derivatives_'];
%
%   To generate the PCA basis (in case you are using a different dataset), use the script
%   testGenerateHarmonicBasis.m, but you need to edit it and change the following
%   variables:
%
% - Resolution of the surfaces: you need to setup the variable params.RES
% (here it is set to 128 x 128). Make sure you have PCA basis of the same
% resolution
%
%% Copyright
%  Author: Hamid Laga (hamid.laga@gmail.com)
%  Last update: 2022/03/23
%%

clf;
close all;
imtool close all;

configure_paths;
currDir = pwd;
dataDir = './sample_surfaces/';
outDir  = './output/';          
fname   = '00000';     % this surface is from DFAUST data input surface to test - it has to be in multiscale representation

%% Loading the surface
load([dataDir fname, '.mat'], 'multiResM');

%% Parameters
params.homeDir        = [pwd '/'];
params.surfaceCodeDir = params.homeDir;
params.currDir        = currDir;

% the resolutions
params.RES = [128, 128];  % For PCA basis only one resolution

params.nres = 1;        % We only need one resolution when using PCA basis

% The parameters for SRNF inversion
params.mysmall  = 0.01;
params.d        = 3;
params.myInner  = @innerS2;
params.cutoff   = 1e-6;
params.stepsize = .1;       % 0.01;
params.itermax  =  5000;    % 10000;     % Max No. of iterations

% Note that initMode is ignored when basis_type = 2 (i.e., PCA basis)
params.initMode = 1;        % 0 - Initialize the inversion using a sphere
                            % 1 - Initialize the inversion using a point on
                            %     the linear path
                            % 2 - Initialize the inversion using the
                            %     previously estimated shape
                            % 3 - starts with a given path (see
                         
params.basis_type = 2;      % 1 - Spherical Harmonic basis (default)
                            % 2 - PCA basis (human shapes). In that case you will need to
                            % go to multiresSRNInversion.m and set the path
                            % to the right PCA basis, and eventulally
                            % adjust the paths to Sebastian's/Ian's codes
                            
params.npca_basis = 20;
params.PCA_BASIS_FNAME_PREFIX = [params.homeDir  'PCABasis/DFAUST_PCA_with_derivatives_'];

params.toVisualize = 0;

%% Normalisation (this is the most time consuming)
multiResM = normalize_translation_multiresSurf(multiResM);
multiResM = normalize_scale_multiresSurf(multiResM);

%% Only use nres resolutions - but make sure that teh nres-th res is the last one
multiResM(params.nres, :,:,:) = multiResM(end, :,:,:);

M_multires = cell(1, params.nres);
for i=1:params.nres
    res = params.RES(i, 1:2);    
    M_multires{i} = resampleSphericalGrid(squeeze(multiResM(i, :,:,:)), res(1), 0);       
end

%% Compute the Qmaps at each resolution
for i=1:params.nres 
    Q_multires{i} = surface2srnf(M_multires{i}); %, Theta);
end

%% Inversion
f_reconstructed_multires = invertSRNF_PCABasis(Q_multires, params);
f_reco = f_reconstructed_multires{end};     % 

%% visualization (Let's just visualize the last resolution)
close all;

i = params.nres;
res= params.RES(i, :);

% Visualize the original surface as projected on the PCA basis
load( [params.PCA_BASIS_FNAME_PREFIX num2str(res(1))  '.mat'], 'B', 'Mu');  
B = B(:, :, 1: min(params.npca_basis, size(B, 3)) );
[cn, fn] = reconstruct_surfaceR3_PCA(mySurf(M_multires{end}), B, res, Mu);

f0 = reshape(fn, [3, res]);
f_hb  = permute(f0,[2,3,1]);
clear('Mu');

h3 = figure(100); clf; 
hsub = subplot(1, 4, 1);
subpos = get(hsub,'position'); 
subpos(1:2) = subpos(1:2) + 75; 
subpos(3:4) = 100; 
uicontrol ('Style','text', 'String', 'Original Surface, with PCA basis','FontSize', 8, 'Position', subpos);
dispSurface(f_hb); 

% The original Surface
hsub = subplot(1, 4, 2);
subpos = get(hsub,'position'); 
subpos(1:2) = subpos(1:2) + 75; subpos(1) = subpos(1) + 100;
subpos(3:4) = 100; 
uicontrol ('Style','text', 'String', 'Original Surface', ...
           'FontSize', 8, 'Position', subpos);
dispSurface(M_multires{i}); 

% Visualize the reconstructed surface
hsub   = subplot(1, 4, 3);
subpos = get(hsub,'position'); 
subpos(1:2) = subpos(1:2) + 75; subpos(1) = subpos(1) + 220;
subpos(3:4) = 100; 
uicontrol ('Style','text', 'String', 'Reconstructed Surface', ...
           'FontSize', 8, 'Position', subpos);
dispSurface(f_reconstructed_multires{i}); 

%% Visualize the reconstruction error
%  Here the error is between f_hb (i.e., the surface as represented with
%  the harmonic basis) and the reconstructed surface. 
%  If you want to see the error with respect to the exact input surface,
%  replace f_hb by M_multires{1}

f_reco  = f_reconstructed_multires{params.nres};
Err = sum((f_hb - f_reco).^2, 3);
hsub = subplot(1, 4, 4);
subpos = get(hsub,'position'); 
subpos(1:2) = subpos(1:2) + 75; subpos(1) = subpos(1) + 320;
subpos(3:4) = 100; 
uicontrol ('Style','text', 'String', 'Reconstruction Error', ...
           'FontSize', 8, 'Position', subpos);
surface(squeeze(f_hb(:,:,1)), squeeze(f_hb(:,:,2)), squeeze(f_hb(:,:,3)), Err); %, ...
colormap jet;
shading interp
light; lighting phong; 
axis equal off;
cameramenu
colorbar

%% Save the figures if needed
% saveas(h, [outDir '/' fname '_' num2str(res(1)) 'x' num2str(res(2)) '_error.fig']);
% save([outDir '/' fname '_reconstructed.mat'] , 'multiResM', 'cn', 'E', 'dE');



