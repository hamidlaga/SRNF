%
% Uses multires surfaces generated with sph wavelet decompositionto compute
% geodesics
%  (this is the correct one)
%
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

drive = '/Users/hamidlaga/Dropbox/'; % PC (edit it according to your own settings)
%global surfaceCodeDir;
%global homeDir;


currDir = pwd;

dataDir = './sample_surfaces/';
outDir  = './sample_results/'; % geodesics_init_with_previous_shape/';
mkdir(outDir);

fname = '00000'; % 'octopus';        % the surface to test

%% Load the surface
load([dataDir fname, '.mat'], 'multiResM');

%% Parameters
params.homeDir        = [drive 'Home/Applications/3DShapeStatistics/code/4DStatistics/SRNF/'];        
params.surfaceCodeDir = params.homeDir;
params.currDir        = currDir;

% the resolutions
params.RES = [ 8  8; ...
       16 16; ...
       16, 16; ...
       16 16;  ...
       16, 16; ...
       25, 25; ...
%        25 25;  ...
%        50 50; ...
       ];

params.LLS = { [2:8];  ...
        [8:36]; ... 
        [36]; ... 
        [36]; ... 
	    [36]; ... 
        [22,23,36]; ... 
        [36]; 
        };
params.nres =  min(size(params.RES, 1), numel(params.LLS));

% The parameters for SRNF inversion
params.mysmall  = 0.01;
params.d        = 3;
params.myInner  = @innerS2;
params.cutoff   = 1e-6;
params.stepsize = .1;       % 0.01;
params.itermax  =  5000; % 10000;     % Max No. of iterations

params.initMode = 1;        % 0 - Initialize the inversion using a sphere
                            % 1 - Initialize the inversion using a point on
                            %     the linear path
                            % 2 - Initialize the inversion using the
                            %     previously estimated shape
                            % 3 - starts with a given path (see
% params.initPath = [outDir 'cat0_cat9_pca/cat0_cat9_geodesic.mat'];

                            
params.basis_type = 1;      % 1 - Spherical Harmonic basis (default)
                            % 2 - PCA basis (human shapes). In that case you will need to
                            % go to multiresSRNInversion.m and set the path
                            % to the right PCA basis, and eventulally
                            % adjust the paths to Sebastian's/Ian's cods
                            
% params.npca_basis = 20;
% params.PCA_BASIS_FNAME_PREFIX = 'PCABasis_HumanShapes/PCABasisHuman_Res'; % 'PCABasis_cat/PCABasisCat_Res';
params.harmonic_basis_homedir = '/Users/hamidlaga/Dropbox/Home/Applications/3DShapeStatistics/Release/InverseQmap/HarmonicBasis/';                            
params.toSave = 1;      % set it to 1 if you want to save the computed gedoesics and figures

%%
multiResM(params.nres, :,:,:) = multiResM(end, :,:,:);
multiResM = multiResM(1:params.nres, :,:,:);
M_multires= prepareMultiResSurface(multiResM, params);

% 
% multiResM = normalize_translation_multiresSurf(multiResM);
% multiResM = normalize_scale_multiresSurf(multiResM);
% 
% [~, O] = alignSurfaceGridWithSphere(squeeze(multiResM(end, :,:,:)), params.surfaceCodeDir, params.currDir); 
% for i=1:size(multiResM, 1),
%     multiResM(i, :,:,:) = rotateGrid(squeeze(multiResM(i, :,:,:)), O);
% end
% 
% M = squeeze(multiResM(end, :,:,:));  % the finest surface
% 
% %% Multi-resolution representation of the surface M_orig
% multiResM(params.nres, :,:,:) = multiResM(end, :,:,:);
% M_multires = cell(1, params.nres);
% for i=1:params.nres,
%     res = params.RES(i, 1:2);    
%     M_multires{i} = resampleSphericalGrid(squeeze(multiResM(i, :,:,:)), res(1), 0);       
% end

%% Compute the Qmaps at each resolution
for i=1:params.nres,
    Q_multires{i} = surface2srnf(M_multires{i});  
end
disp('Qmap computed, press anykey to continue ...');

%% Now processing
f_reconstructed_multires = invertSRNF(Q_multires, params);

%% The surface as represented by the spherical harmonic basis
load( [params.harmonic_basis_homedir  'harmonic_basis_Res' num2str(res(1)) '_L' num2str(params.LLS{end}(end)) '.mat']); 
[c_new, f_new] = reconstruct_surfaceR3( mySurf(M_multires{end}), @innerS2, B, res);
f_hb = permute(reshape( f_new, [3,res(1), res(2)]), [2, 3, 1]);

%% Save the reconstructed surface (if needed)
% save([outDir '/' fname '_reconstructed.mat'] , 'multiResM', 'cn', 'E', 'dE');

%% visualization
for i=params.nres:params.nres,
    res= params.RES(i, :);
        
    h1 = figure(100); clf; 
    dispSurface(M_multires{i}); 
   
    h2 = figure(200); clf;
    dispSurface(f_reconstructed_multires{i}); 
    
    disp('Press anykey to continue ...');
    pause; % (.1);
end
h3 = figure(300); clf;
dispSurface(f_hb); 


%% visualize the reconstruction error (if needed)
f_reco  = f_reconstructed_multires{params.nres};
Err = sum((f_hb - f_reco).^2, 3);
h  = figure(400); clf;

surface(squeeze(f_hb(:,:,1)), squeeze(f_hb(:,:,2)), squeeze(f_hb(:,:,3)), Err); %, ...
colormap jet;
shading interp
light; lighting phong; 
axis equal off;
cameramenu
colorbar
% saveas(h, [outDir '/' fname '_' num2str(res(1)) 'x' num2str(res(2)) '_error.fig']);

