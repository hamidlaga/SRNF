%
% Uses multires surfaces generated with sph wavelet decomposition
%  (this is the correct one)
% This code uses PCA basis -  e.g., for human body shapes
%
% REQUIREMENTS
% - Code in code/Surface
%
% - Harmonic Basis for generic surfaces - generateMultipleHarmonicBasis.m
%   It is recommended that you pre-generate them and just load them.
%
% - PCA basis for human shapes
%
% - Compile the mex files
%       mex innerS2.cpp myVector.cpp
%       mex  map2qnew.cpp myVector.cpp   
%       mex  gradInvQnew.cpp myVector.cpp
%
% Preparation of the surfaces:
% - Spherical parameterization
% - Multiresolution representation using spherical harmonics
% - 
%

clf;
close all;
imtool close all;

configure_paths;

drive   = '/Users/hamidlaga/Dropbox/'; % edit it according to your own settings
currDir = pwd;
dataDir = './sample_surfaces/';
outDir  = './sample_results/';         % geodesics_init_with_previous_shape/';
mkdir(outDir);

fname = '00000'; %'00111'; % 'michael10';        % the surface to test

%% Loading the surface
load([dataDir fname, '.mat'], 'multiResM');

%% Parameters
params.homeDir = [drive 'Home/Applications/3DShapeStatistics/code/4DStatistics/SRNF/'];
params.surfaceCodeDir = params.homeDir;
% the resolutions
params.RES = [32, 32];  % For PCA basis only one resolution

params.nres = 1; 

% The parameters for SRNF inversion
params.mysmall  = 0.01;
params.d        = 3;
params.myInner  = @innerS2;
params.cutoff   = 1e-6;
params.stepsize = .1;       % 0.01;
params.itermax  =  5000;    % 10000;     % Max No. of iterations

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
params.PCA_BASIS_FNAME_PREFIX = '/Users/hamidlaga/Dropbox/Home/Applications/3DShapeStatistics/code/4DStatistics/PCABasis/DFAUST_PCA_with_derivatives_';

params.toVisualize = 0;

%%
% load('./sample_surfaces/Lthalamus.mat');
% S = squeeze(data(:, :, :, 1));

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

%% Visualize the initial surface
if 1, %params.toVisualize
    M = squeeze(multiResM(end, :,:,:));  % the finest surface
    figure, dispSurfR3(mySurf(M), [size(M, 1), size(M, 2)], 3);
    
    % Visualization (the inverted surface)
    for i=params.nres:params.nres,
        res= params.RES(i, :);

        % The original surface
        h1 = figure(100); clf; 
        dispSurface(M_multires{i}); 

        % The reconstructed surface
        h2 = figure(200); clf;
        dispSurface(f_reconstructed_multires{i}); 

        disp('Press anykey to continue ...');
        pause; % (.1);
    end

    %% The surface as represented by the PCA bases
    % load( [params.harmonic_basis_homedir  'harmonic_basis_Res' num2str(res(1)) '_L' num2str(LLS{end}(end)) '.mat']); 
    % [c_new, f_new] = reconstruct_surfaceR3( mySurf(M_multires{end}), @innerS2, B, res);
    % f_hb = permute(reshape( f_new, [3,res(1), res(2)]), [2, 3, 1]);

    load( [params.PCA_BASIS_FNAME_PREFIX num2str(res(1))  '.mat'], 'B', 'Mu');  
    B = B(:, :, 1: min(params.npca_basis, size(B, 3)) );
    [cn, fn] = reconstruct_surfaceR3_PCA(mySurf(M_multires{end}), B, res, Mu);
    clear('Mu');

    %% Visualize the original surface
    figure(300), 
    dispSurfR3(fn, res, 3); 

    %% visualize the reconstruction error (if needed) as a color map on the original surface
    f0 = reshape(fn, [3, res]);
    f0  = permute(f0,[2,3,1]);
    Err = sum((f0 - f_reco).^2, 3);

    h   = figure(400); clf;
    surface(squeeze(f0(:,:,1)), squeeze(f0(:,:,2)), squeeze(f0(:,:,3)), Err); %, ...
    colormap jet;
    shading interp
    light; lighting phong; 
    axis equal off;
    cameramenu
    colorbar

end

%% Save the figures if needed
% saveas(h, [outDir '/' fname '_' num2str(res(1)) 'x' num2str(res(2)) '_error.fig']);
% save([outDir '/' fname '_reconstructed.mat'] , 'multiResM', 'cn', 'E', 'dE');



