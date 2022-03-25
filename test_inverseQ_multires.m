%% README
% This script demonstrates how to use the SRNF inversion procedure.
%
% In this example, we show how to compute the invert of an SRNF map when an
% an initial surface is not available. In this case, the initialization will
% be a sphere.
%
% The algorithm uses a multiscale representation of surfaces. The examples
% provided here have been generated using sphercial wavelet transform
% 
% REQUIREMENTS
% - Harmonic Basis for generic surfaces 
%   the script for this is testGenerateHarmonicBasis.m, 
%   Its header contains an explanation on how to use it.
%
%   It is recommended that you pre-generate the basis and just load them here.
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
% - Multiscale representation using spherical harmonics
% If you plan to initialize the inversion with a sphere, you need to work with 
% multiscale surfaces. The examples provided here have been generated using 
% spherical wavelet transform. The code will be made available soon. 
%
%% Example
% The script below takes a surface, computes its SRNF and then tries to
% recover the initial surface by SRNF inversion. It uses a sphere as an
% initialization
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
% testGenerateHarmonicBasis.m, but you need to edit it and change the following
% variables:
%        -- res to the resolutions you want and 
%        -- L: to the frequencies you want
%        -- outDir: the output directort where the generated harmonic basis
%        will be saved
% Note that for this example, we are using the following resolutions and frequencues::
%         8 x 8: LL from 2 to 8
%         16 x 16: LL from 8 to 36
%         25 x 25: LL of 22, 23 and 36
% You may need more if dealing with high resolution surfaces. 
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
fname   = '00000';    % the input surface to test - it has to be in multiscale representation

%% Load the surface
load([dataDir fname, '.mat'], 'multiResM');

%% Parameters
params.homeDir        = [pwd '/'];
params.surfaceCodeDir = params.homeDir;
params.currDir        = currDir;

% the resolutions (you can increase further the resolution, but will be
% slower)
params.RES = [ 8  8; ...
       16 16; ...
       16, 16; ...
       16 16;  ...
       16, 16; ...
       25, 25; ...
%        25 25;  ...
%        50 50; ...
       ];

params.LLS = { [2:8];  ... % each row should correspond to one row in params.RES above
        [8:36]; ... 
        [36]; ... 
        [36]; ... 
	    [36]; ... 
        [22,23,36]; ... 
%        [36]; ...
%        [36]; ...
        };
params.nres =  min(size(params.RES, 1), numel(params.LLS));

% The parameters for SRNF inversion
params.mysmall  = 0.01;
params.d        = 3;
params.myInner  = @innerS2;
params.cutoff   = 1e-6;
params.stepsize = .1;       % 0.01;
params.itermax  =  5000;    % Max No. of iterations

params.initMode = 1;        % 0 - Initialize the inversion using a sphere
                            % 1 - Initialize the inversion using a point on
                            %     the linear path
                            % 2 - Initialize the inversion using the
                            %     previously estimated shape
                            % 3 - starts with a given path (see
                            
params.basis_type = 1;      % 1 - Spherical Harmonic basis (default)
                            
params.harmonic_basis_homedir = [params.homeDir '../HarmonicBasis/'];  
params.toSave = 1;      % set it to 1 if you want to save the computed gedoesics and figures (not used)

params.cn = [];         % empty for a sphere - otherwise, specify the harmonic rep of the initial surface

%%
multiResM(params.nres, :,:,:) = multiResM(end, :,:,:);
multiResM  = multiResM(1:params.nres, :,:,:);
M_multires = prepareMultiResSurface(multiResM, params);

%% Compute the Qmaps at each resolution
for i=1:params.nres,
    Q_multires{i} = surface2srnf(M_multires{i});  
end
disp('Qmap computed ...');

%% Now processing
f_reconstructed_multires = invertSRNF(Q_multires, params);

%% The surface as represented by the spherical harmonic basis
res = params.RES(end, :);
load( [params.harmonic_basis_homedir  'harmonic_basis_Res' num2str(res(1)) '_L' num2str(params.LLS{end}(end)) '.mat']); 
[c_new, f_new] = reconstruct_surfaceR3( mySurf(M_multires{end}), @innerS2, B, res);
f_hb = permute(reshape( f_new, [3,res(1), res(2)]), [2, 3, 1]);

%% Save the reconstructed surface (if needed)
% if params.toSave
%     disp([ 'Reconstructed surface is save into ' outDir '/' fname '_reconstructed.mat']);
%     save([outDir '/' fname '_reconstructed.mat'] , 'multiResM', 'cn', 'E', 'dE');
% end

%% visualization (Let's just visualize the last resolution)
close all;
h3 = figure(100); clf; %title('Original Surface, with harmonic bases');
hsub = subplot(1, 4, 1);
subpos = get(hsub,'position'); 
subpos(1:2) = subpos(1:2) + 75; 
subpos(3:4) = 100; 
uicontrol ('Style','text', 'String', 'Original Surface, with harmonic basis','FontSize', 8, 'Position', subpos);
dispSurface(f_hb); 

i  = params.nres;
res= params.RES(i, :);

% Original Surface
hsub = subplot(1, 4, 2);
subpos = get(hsub,'position'); 
subpos(1:2) = subpos(1:2) + 75; subpos(1) = subpos(1) + 100;
subpos(3:4) = 100; 
uicontrol ('Style','text', 'String', 'Original Surface', ...
           'FontSize', 8, 'Position', subpos);
dispSurface(M_multires{i}); 

% Reconstructed surface
hsub = subplot(1, 4, 3);
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


