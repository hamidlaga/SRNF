%% README
%% Geodesic shooting for the human shapes
% This script demonstrates how to compute geodesics between human body shapes using SRNF inversion procedure.
% This example uses PCA basis computed on a subset of the [DFAUST](https://dfaust.is.tue.mpg.de/) dataset.
%
% It can also be used for deformation transfer using geodesic shooting (see the detailed description below).  
%
% This code uses spherical grids of 128 x 128.
%
% REQUIREMENTS
% - PCA basis of the class of shapes. This example uses PCA basis computed on a subset of the [DFAUST](https://dfaust.is.tue.mpg.de/) dataset.  
%   If you are using a differen human dataset (or class of shapes) you need
%   to generate its corresponding PCA bases. Use the script
%   testGeneratePCABasis.m
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
% - Elastic registration: this script assumes that the surfaces are
% pre-registered. Please refer to section "Elastic registration and geodesics" in https://sites.google.com/view/hamidlaga/3d4d-shape-analysis 
%
%% Example
%% Computing the geodesic between two surfaces S1 and S2
% fname1 corresponds to S1
% fname2 corresponds to S2
% fanme3 corresponds to S1
% 
% The code will then should a geodesic from S1 to S2 and transfers it to S1
%
%% Transfer the deformation (S1, S2) to another surface S3, which is on the same pose as S1
% fname1 corresponds to S1
% fname2 corresponds to S2
% fanme3 corresponds to S3
%
% The code will then should a geodesic from S1 to S2 and transfers it to S3
%
% Note that in the example below you may need to edit the paths depending
% on where you saved the code and the data.
%
% Requirements
% - PCA basis should be pre-generated. These should be generated and saved in the folder
% ./PCABasis
%   If you want to use a different folder, please change the variable
%     params.harmonic_basis_homedir = [pwd '/PCABasis/'];
%
% These PCA basis should correspond to the class of shapes you are dealing
% with. This example uses PCA basis computed on a subset of the [DFAUST](https://dfaust.is.tue.mpg.de/) dataset.
% 
% To generate the PCA basis for this example, use the script
% testGeneratePCA.m, but you need to edit it some of the lines of the code. 
%
%% Copyright
%  Author: Hamid Laga (hamid.laga@gmail.com)
%  Last update: 2022/03/23
%%

clear all;
close all;
configure_paths;

mysmall = 0.01;

global surfaceCodeDir;
currDir        = pwd;
homeDir        = [pwd '/']; 
surfaceCodeDir = [homeDir '../Surface/'];

dataDir = './sample_surfaces/';
outDir  = './output/';

fname1 = '50002_chicken_wings_00000'; 
fname2 = '50002_chicken_wings_00077';
fname3 = fname1;  

% The desired resolution
res = 128;

%% Parameters
params.homeDir        = homeDir;
params.surfaceCodeDir = surfaceCodeDir;
params.currDir        = currDir;

% the resolutions (you can increase further the resolution, but will be
% slower)
params.RES    = [res, res]; 
params.LLS{1} = [25];
params.nres   =  min(size(params.RES, 1), numel(params.LLS));

% The parameters for SRNF inversion
params.mysmall  = 0.0;
params.d        = 3;
params.myInner  = @innerS2;
params.cutoff   = 1e-6;
params.stepsize = .1; 
params.itermax  =  10000;   % Max No. of iterations 

params.initMode = 1;        % 0 - Initialize the inversion using a sphere
                            % 1 - Initialize the inversion using a point on
                            %     the linear path
                            % 2 - Initialize the inversion using the
                            %     previously estimated shape
                            % 3 - starts with a given path (see
                            
params.basis_type = 2;      % 1 - Spherical Harmonic basis (default)
                            % 2 for PCA,
                            
params.npca_basis = 100;    % No. of PCA basis
params.PCA_BASIS_FNAME_PREFIX = './PCABasis/DFAUST_PCA_with_derivatives_';
 
params.toSave = 1;      % set it to 1 if you want to save the computed gedoesics and figures (not used)

params.cn = [];         % empty for a sphere - otherwise, specify the harmonic rep of the initial surface

params.nsteps = 5;      % Number of samples along the geodesic

%% Load two registered surfaces (fname1 and fname2)
load([dataDir fname1, '.mat'], 'M' );  
M(:, end,:)= [];
M(:, 1, :) = [];
M(end, :,:)= [];
M(1, :,:)  = [];
multiResM1(1, :, :, :) = resampleSphericalGrid(M, res(1), 0);

load([dataDir fname2, '.mat'], 'M' );  
M(:, end,:)= [];
M(:, 1,  :)= [];
M(end, :,:)= [];
M(1, :,  :)= [];
multiResM2(1, :,:,:) = resampleSphericalGrid(M, res(1), 0);

%% Load the shape to deform
multiResM3 = multiResM1;

% normalize for translation
multiResM1 = normalize_translation_multiresSurf(multiResM1);
multiResM2 = normalize_translation_multiresSurf(multiResM2);
multiResM3 = normalize_translation_multiresSurf(multiResM3);

% normalize for scale
multiResM1 = normalize_scale_multiresSurf(multiResM1);
multiResM2 = normalize_scale_multiresSurf(multiResM2);
multiResM3 = normalize_scale_multiresSurf(multiResM3);

%% Rigid Alignment of M2 on M1
M1 = squeeze(multiResM1(end, :,:,:));        
M2 = squeeze(multiResM2(end, :,:,:));        

[~, O] = rigidAlign(M1, M2, surfaceCodeDir, currDir);
for i=1:size(multiResM2, 1),
    multiResM2(i, :,:,:) = rotate3D(squeeze(multiResM2(i, :,:,:)), O);
end

% Rigid alignment of M3 on M1
M3 = squeeze(multiResM3(end, :,:,:));       
[~, O] = rigidAlign(M1, M3, surfaceCodeDir, currDir);
for i=1:size(multiResM3, 1),
    multiResM3(i, :,:,:) = rotate3D(squeeze(multiResM3(i, :,:,:)), O);
end

M1_multires = cell(1, params.nres);
M2_multires = cell(1, params.nres);

for i=1:params.nres,
    res = params.RES(i, 1:2);      
    K = resampleSphericalGrid(squeeze(multiResM1(i, :,:,:)), res(1), 0);
    [A1,~,~] = area_surf_closed(K);   A1 = 1; 
    M1_multires{i} = K / sqrt(A1); 
    
    K = resampleSphericalGrid(squeeze(multiResM2(i, :,:,:)), res(1), 0);
    [A2,~,~] = area_surf_closed(K); A1 = 1;
    M2_multires{i} = K / sqrt(A2); 
    
    K = resampleSphericalGrid(squeeze(multiResM3(i, :,:,:)), res(1), 0);
    [A3,~,~] = area_surf_closed(K); A1 = 1;
    M3_multires{i} = K / sqrt(A3); 
    
end

%% Compute the Qmaps at each resolution
Q1_multires = cell(1, params.nres);
Q2_multires = cell(1, params.nres);
Q3_multires = cell(1, params.nres);
for i=1:params.nres,
    res = params.RES(i, 1:2);    
    [Theta,Phi] = genGridSphr(res, params.mysmall);    
    sinTheta  = genSinPhi(Theta, params.mysmall);
    N       = prod(res);
    d       = 3;

    res = params.RES(i, :);
    mySmall = 0.01;
    Q11 = surface2qnew_qie(mySurf(M1_multires{i}, 0), res, mySmall); 
    Q1_multires{i} = permute(reshape(Q11, [3, res(1), res(2)]), [2 3 1]);
    
    Q22 = surface2qnew_qie(mySurf(M2_multires{i}, 0), res, mySmall); 
    Q2_multires{i} = permute(reshape(Q22, [3, res(1), res(2)]), [2 3 1]);
    
    Q33 = surface2qnew_qie(mySurf(M3_multires{i}, 0), res, mySmall); 
    Q3_multires{i} = permute(reshape(Q33, [3, res(1), res(2)]), [2 3 1]);
end

%% Shooting vectors
%  shooting vector between M1 and M2
V_multiRes = cell(1, params.nres);
for i = 1:params.nres,   
    f1   = squeeze(M1_multires{i});
    f2   = squeeze(M2_multires{i});
    V_multiRes{i} = f2 - f1;
end

% Shooting vector from Q1 to Q2
QV_multiRes = cell(1, params.nres);
for i = 1:params.nres,   
    q1   = squeeze(Q1_multires{i});
    q2   = squeeze(Q2_multires{i});
    QV_multiRes{i} = q2 - q1;
end

%% Shooting in the geodesic direction
S      = linspace(0, 1, params.nsteps); % [0 .25 .5 .75, 1 , 1.2, 1.5, 1.75, 2]; % [0 .25, .5, .75 , 1.0, 1.25, 1.5, 1.75, 2];  % I am going to test the extrapolations
nsteps = numel(S);
Qpath_multires = cell(1, params.nsteps);
Fpath_multires = cell(1, params.nsteps);

for j=1:nsteps
    Q = cell(1, params.nres);
    F = cell(1, params.nres);
    s = S(j);
        
    for i = 1:params.nres,   
        q    = squeeze(Q3_multires{i});
        q1   = squeeze(Q1_multires{i});
        q2   = squeeze(Q2_multires{i});
        qV   = QV_multiRes{i};      
        Q{i} = transportSRNFdeformation(q1, q2, s, q); % q + s * qV; % 
     
        f    = squeeze(M3_multires{i});
        fV   = V_multiRes{i};
        F{i} = f + s * fV;
    end
    
    Qpath_multires{j} = Q;
    Fpath_multires{j} = F;
    
end

%% the inversion of the entire path
geodesicPath{1}   = M3_multires;
res = params.RES(end, :);
L   = params.LLS{end};
l   = L(end);
    
cn_multires = [];
Elinpath  = zeros(nsteps,1);        % energy of the linear path
Egeodpath = zeros(nsteps,1);        % energy of the geodesic path

geodesicPath{params.nsteps} = M2_multires; 
for ii = 2: params.nsteps, 
    % this is the Q to invert
    Q_multires  = Qpath_multires{ii};
    
    %% inversion (comment / uncomment to select the appropriate option)      
    % starting with the linear path
    f0_multires = Fpath_multires{ii}; % linearPath{ii}; % geodesicPath{ii-1}; %    
    f0 = f0_multires{1};   
    
    res = params.RES(1, :);
    
    if params.basis_type == 1, % harmonic basis (not used in this example)
        load( ['HarmonicBasis/harmonic_basis_Res' num2str(res(1)) '_L' num2str(l) '.mat'], 'B');  
        [cn, f0n] = reconstruct_surfaceR3(f0, params.myInner, B, res);
    else % PCA basis
        load( [params.PCA_BASIS_FNAME_PREFIX  num2str(res(1)) '.mat'], 'B', 'Mu');  
        B = B(:,:, 1: min(params.npca_basis,  size(B, 3) ));
        [cn, f0n] = reconstruct_surfaceR3_PCA(mySurf(f0), B, res, Mu);
        clear('Mu');
    end
    clear('B');
  
    Qs = Q_multires;
    LL = params.LLS;
    RESn = params.RES;
    tic
    [f_multires, cn_multires, E, dE] = multiresSRNFInversion(Q_multires,  cn, params); %% params.RES(1:params.nres, :), LLS,
    toc;
    Elinpath(ii)  = E{1}(1);
    Egeodpath(ii) = E{end}(end);
    
    disp([Elinpath(ii), Egeodpath(ii)]);
    
    % Save
    geodesicPath{ii} = f_multires;
    
    % disp('Press anykey to continue ...');
    % pause(.1);
    
end

%% Visualizing the computed geodesics
% Visualizing the geodesic at hires
Fgeod = [];
R0 = getRotationMatrix([1, 0, 0], 65,  0);

for i=1: nsteps,
    R  = normalize_scale_multiresSurf(geodesicPath{i});    
    Fgeod(:,:,i,:) = rotate3D(rotateGrid(R{end}, inv(O)), R0);
end

fig_geodesics = 2000;
hgeod = figure(fig_geodesics); clf;
DisplayGeod(Fgeod, fig_geodesics,30,30);

% Visualizing the linear path at hires
nres = 1;
O  = eye(3, 3);
R0 = eye(3, 3);
M1 = squeeze(multiResM1(1, :,:,:)); 
M2 = squeeze(multiResM2(1, :,:,:)); 
Flin = [];                              % 
for i=1:nsteps,
    s = (i-1)/(nsteps-1);
    Flin(:,:,i,:) = rotate3D(rotateGrid(M1*(1-s)+ M2*s, inv(O)), R0);    
end

fig_linear = 4000;
hlin = figure(fig_linear); clf;
DisplayGeod(Flin, fig_linear,30,30);

%% Uncomment if you want to save the results
% F = Fgeod;
% save([outDir fname '_geodesic.mat'], 'F');
% saveas(hgeod, [outDir fname '_geodesic.fig']);
% print([outDir fname '_geodesic.eps'],'-depsc');

% F = Flin;
% save([outDir fname '_linearpath.mat'], 'F');
% saveas(hlin, [outDir fname '_linearpath.fig']);
% print([outDir fname '_linearpath.eps'],'-depsc');
    
% fig_energy = 5000;
% h3 = figure(fig_energy); clf;
% hold on;
% plot(Egeodpath, '-xr', 'LineWidth', 3);
% plot(Elinpath, '-ob', 'LineWidth', 3);
% set(gca, 'FontSize', 32);
% saveas(h3, [outDir fname '_pathenergies.fig']);
% print([outDir fname '_pathenergies.eps'],'-depsc');
% 
% 
