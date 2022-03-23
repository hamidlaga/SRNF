%% README
% This file demonstrate how to use the SRNF inversion procedure.
%
% In this example, we show how to compute the invert of an SRNF map when an
% initial surface is give. Importantly, the initial surface should be close
% to the surface we want to recover.
% 
% REQUIREMENTS
% - Harmonic Basis for generic surfaces 
%   the script for this is testGenerateHarmonicBasis.m, 
%   Its header contains an explanation on how to use it.
%
%   It is recommended that you pre-generate the basis and just load them here.
%
% - PCA basis for human shapes
%   If you are dealing with a specific class of shapes (e.g., humans), it is better to use
%   a basis that is specific to that class, e.g., PCA basis computed on some training data.  
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
%% Example
% The script below takes a surface, computes its SRNF and then tries to
% recover the initial surface by SRNF inversion. It uses the linear mean in the
% space of surfaces to initialize the inversions. It assumes that the mean
% is close to the surface thus the code does not require a muttiresolution representation
% of the surfaces to perform the inversion. 
%
% Note that in the example below you may need to edit the paths depending
% on where you saved the code and the data.
%
% Requirements
% - harmonic basis should be pre-generated. These should be generated and saved in the folder
% ./HarmonicBasis
%   If you want to use a different folder, please change the variable
%     params.harmonic_basis_homedir = [pwd '/HarmonicBasis/'];
%
% To generate the harmonic basis for this example, use the script
% testGenerateHarmonicBasis.m, but you need to it and change the following
% variables:
%        -- res to 50 since this example uses spherical surfaces of
%        resolution 50 x 50
%        -- L:the default one is set from 2 to 36. However, for this
%        example we only need L = [22,23,36]
%        -- outDir: the output directort where the generated harminoc basis
%        will be save
%
%%
% Copyright
% Author: Hamid Laga (hamid.laga@gmail.com)
% Latest update: 2022/03/22
%
clf;
close all;
imtool close all;

configure_paths;
currDir = pwd;

dataDir   = './sample_surfaces/';
outputDir = './output/';                 % output directory

%% Load the mat file that contains the surface to invert and the initial surface (for initialization)
load([dataDir './surfaces.mat'], 'F', 'muF');
nSurfaces = 1; 

% The resolution of the surfaces
res(1) = size(muF, 1);
res(2) = size(muF, 2);

% Let's work with only one resolution (since we initialise the inversion
% with the linear mean
nres      = 1;                              % Just one resolution
multiResM = zeros(nres, res(1), res(2), 3); % This will hold the multires surface (just one resolution here)
multiResM(1, :,:,:) = resampleSphericalGrid(F, res(1), 0);

%% Parameters
params.homeDir        = [pwd '/']; 
params.surfaceCodeDir = params.homeDir;
params.currDir        = currDir;

params.RES = res;   % the resolutions
params.LLS = {      % The levels of the harmonic decomposition we will consider
        [22,23,36]
        };
params.nres = min(size(params.RES, 1), numel(params.LLS));

L   = params.LLS{end};
l   = L(end);

% The parameters for SRNF inversion
params.mysmall  = 0.01;
params.d        = 3;
params.myInner  = @innerS2;
params.cutoff   = 1e-6;
params.stepsize = .1;       % 0.01;
params.itermax  =  5000;    % 10000;     % Max No. of iterations

% How to initialize the inversion
params.initMode = 2;        % 0 - Initialize the inversion using a sphere
                            % 1 - Initialize the inversion using a point on
                            %     the linear path
                            % 2 - Initialize the inversion using the
                            %     previously estimated shape (the provided
                            %     linear mean in this case)
                            %     
                            % 3 - starts with a given path. If you choose this, make sure you provide the path
                            %     e.g., params.initPath = [outDir 'cat0_cat9_pca/cat0_cat9_geodesic.mat'];
                            
params.basis_type = 1;      % 1 - Spherical Harmonic basis (default)
                            % 2 - PCA basis (human shapes). In that case you will need to
                            % go to multiresSRNInversion.m and set the path
                            % to the right PCA basis, and eventulally
                            % adjust the paths to Sebastian's/Ian's cods

% Uncomment the following two if  you are using PCA basis                            
% params.npca_basis = 20;   % no. of PCA basis
% params.PCA_BASIS_FNAME_PREFIX = 'PCABasis_HumanShapes/PCABasisHuman_Res'; % Path the PCA basis where they are stored
% 'PCABasis_cat/PCABasisCat_Res';

% Path to where the pre-generated harmonic basis are stored
params.harmonic_basis_homedir = [pwd '/HarmonicBasis/'];                           
params.toSave = 0;      % set it to 1 if you want to save the computed gedoesics and figures

% Initial surface
params.cn = [];     % empty for a sphere - otherwise, specify the harmonic rep of the initial surface

%% Initial surface - projecting it onto the harmonic basis
res = params.RES;
initialM(1, :,:,:) = muF;
intialM_multires= prepareMultiResSurface(initialM, params);

load( [params.harmonic_basis_homedir 'harmonic_basis_Res' num2str(res(1)) '_L' num2str(l) '.mat'], 'B');  
[cn0, f0n] = reconstruct_surfaceR3(mySurf(intialM_multires{1}, 0), params.myInner, B, res);
dispSurfR3(f0n, res, 3);

params.cn = cn0;    %the harmonic coeffs of the initial surface

%% Compute the Qmaps of the input surfaces
M_multires= prepareMultiResSurface(multiResM, params);
for i=1:nres,
    res     = params.RES(i, :);
    mySmall = 0.01;
    Q11     = surface2qnew_qie(mySurf(M_multires{i}, 0), res, mySmall); % , sinTheta, d, N, res);    
    Q_multires{i} = permute(reshape(Q11, [3, res(1), res(2)]), [2 3 1]);
end

%% Now let's try to invert the SRNF - the result is in M
f_reconstructed_multires = invertSRNF(Q_multires, params);
M = f_reconstructed_multires{1};

% Save the reconstructed surface 
save([outputDir num2str(1) '.mat'], 'M');

%% Let's measure the reconstruction errors
% The surface as represented by the spherical harmonic basis
[c_new, f_new] = reconstruct_surfaceR3( mySurf(M_multires{end}), @innerS2, B, res);
f_hb = permute(reshape( f_new, [3,res(1), res(2)]), [2, 3, 1]);

%% Save the reconstructed surface (uncomment the following line if needed)
% save([outDir '/' fname '_reconstructed.mat'] , 'multiResM', 'cn', 'E', 'dE');

%% visualization (this could be significantly improved)
% The surface as represented with the harmonic basis
close all;
h3 = figure(100); clf; %title('Original Surface, with harmonic bases');
hsub = subplot(1, 4, 1);
subpos = get(hsub,'position'); 
subpos(1:2) = subpos(1:2) + 75; 
subpos(3:4) = 100; 
uicontrol ('Style','text', 'String', 'Original Surface, with harmonic basis','FontSize', 8, 'Position', subpos);
dispSurface(f_hb); 

res= params.RES(1, :);

% Original Surface
hsub = subplot(1, 4, 2);
subpos = get(hsub,'position'); 
subpos(1:2) = subpos(1:2) + 75; subpos(1) = subpos(1) + 100;
subpos(3:4) = 100; 
uicontrol ('Style','text', 'String', 'Original Surface', ...
           'FontSize', 8, 'Position', subpos);
dispSurface(M_multires{1}); 

% Reconstructed surface
hsub = subplot(1, 4, 3);
subpos = get(hsub,'position'); 
subpos(1:2) = subpos(1:2) + 75; subpos(1) = subpos(1) + 220;
subpos(3:4) = 100; 
uicontrol ('Style','text', 'String', 'Reconstructed Surface', ...
           'FontSize', 8, 'Position', subpos);
dispSurface(f_reconstructed_multires{1}); 

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
% saveas(h, [outDir '/' fname '_' num2str(res(1)) 'x' num2str(res(2)) '_error.fig']);

