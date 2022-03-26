%% tests geodesic shooting
% Can be used for computing geodesics between 2 registered surface.
% It can also be used to transfer deformation, say (S1 --> S2) will be
% transferred to S3 which is inthe same pose as S1
%
% Hamid Laga
% Last update: 2022/3/26
%

clear all;
close all;
imtool close all;

params.mysmall = 0.01;

params.homeDir        = [pwd '/']; 
params.surfaceCodeDir = params.homeDir; % [params.homeDir '../Surface/'];
params.currDir        = pwd;

dataDir = './sample_surfaces/ForGeodesicShooting/';
outDir  = './output/'; 

params.RES = [25, 25];  % the resolution
%     ...
%         25, 25; ...
%         50, 50; ...
%         50, 50; ...
%         50, 50]; 

params.LLS{1} = [36];  % The scales (frequencies)
% LLS{1} = [2:8]; 
% LLS{2} = [8:16];
% LLS{3} = [16:2:22];
% LLS{4} = [22]; 
% LLS{5} = [22];
% LLS{6} = [36]; 

params.nres = min([numel(params.LLS), size(params.RES, 1)]);

params.d        = 3;
params.myInner  = @innerS2;
params.cutoff   = 1e-6;
params.stepsize = .1;       %0.01;
params.itermax  = 5000;     % Around 200 should be sufficient
params.basis_type = 1;      % spherical harmonics

params.harmonic_basis_homedir = [params.homeDir '../HarmonicBasis/'];

% Steps along the geodesic - for extrapolation, put the range outside [0, 1]
S = linspace(0, 1, 5);      % Example of \eExtrapolation beyong the last shape S = [0 .25 .5 .75 1.0 1.25];  

fname1  = 'Cat1_1_Cat1_4';    % 'Cat1_1_Cat1_2';   % 'Cat1_1_horse-gallop0';    % This file contains the source and target shape
fname2  = 'Cat1_1_Cat1_4';    % 'Cat1_1_Cat1_2';   % 'Cat1_1_horse-gallop0'; %    % For defomration transfer, specify another file
outfname = [fname1 '_' fname2]; %  horse-gallop0_deformed_to_Cat1_2';

%% Load the registered surfaces
load([dataDir fname1, '.mat'], 'multiResM1', 'multiResM2');  
% For deformation transfer, load the 3rd 3D model (multiResM3)- it shoudl be in the
% same pose as multiResM1. The deformation betwee multiResM1 and multiResM2
% will be transferred to multiResM3
%
% For geodesics, just set multiResM3 = MM.multiResM1
%
% MM = load([dataDir fname2, '.mat']);  
multiResM3 = multiResM1; 

% normalize for translation
multiResM1 = normalize_translation_multiresSurf(multiResM1);
multiResM2 = normalize_translation_multiresSurf(multiResM2);
multiResM3 = normalize_translation_multiresSurf(multiResM3);

% normalize for scale
multiResM1 = normalize_scale_multiresSurf(multiResM1);
multiResM2 = normalize_scale_multiresSurf(multiResM2);
multiResM3 = normalize_scale_multiresSurf(multiResM3);

%% remove the poles from narrow features (this is done manually)
%  use the function getRandomRotationMatrix to get a random rotation matrix
% O = eye(3, 3);
% O = [0.5794   -0.5941   -0.5580; ...
%     0.6410   -0.0907    0.7622; ...
%    -0.5034   -0.7993    0.3283];
% O = getRandomRotationMatrix;

% M2 = squeeze(multiResM2(end, :,:,:)); 
% [~, O] = alignSurfaceGridWithSphere(M2, surfaceCodeDir, currDir); 
% for i=1:size(multiResM1, 1),
%     F = squeeze(multiResM1(i, :,:,:));
%     multiResM1(i, :,:,:) = rotateGrid(F, O);
% end
% for i=1:size(multiResM2, 1),
%     F = squeeze(multiResM2(i, :,:,:));
%     multiResM2(i, :,:,:) = rotateGrid(F, O);
% end
% for i=1:size(multiResM3, 1),
%     F = squeeze(multiResM3(i, :,:,:));
%     multiResM3(i, :,:,:) = rotateGrid(F, O);
% end

%% [Uncomment if you want to do Rigid Alignment] Rigid alignment of M2 on M1
% % M2 on M1
% M1 = squeeze(multiResM1(end, :,:,:)); 
% M2 = squeeze(multiResM2(end, :,:,:));        
% [~, O] = rigidAlign(M1, M2, surfaceCodeDir, currDir);
% for i=1:size(multiResM2, 1),
%     multiResM2(i, :,:,:) = rotate3D(squeeze(multiResM2(i, :,:,:)), O);
% end
% 
% % M3 on M1
% M3 = squeeze(multiResM3(end, :,:,:));        
% [~, O] = rigidAlign(M1, M3, surfaceCodeDir, currDir);
% for i=1:size(multiResM3, 1),
%     multiResM3(i, :,:,:) = rotate3D(squeeze(multiResM3(i, :,:,:)), O);
% end

%% Surfaces at the finest resolutions     
M1 = squeeze(multiResM1(end, :,:,:)); 
M2 = squeeze(multiResM2(end, :,:,:));  
M3 = squeeze(multiResM3(end, :,:,:)); 

%% Multi-resolution representation of the surface M_orig
multiResM1(params.nres, :,:,:) = multiResM1(end, :,:,:);
multiResM2(params.nres, :,:,:) = multiResM2(end, :,:,:);
multiResM3(params.nres, :,:,:) = multiResM3(end, :,:,:);
M1_multires = cell(1, params.nres);
M2_multires = cell(1, params.nres);
M3_multires = cell(1, params.nres);

for i=1:params.nres,
    res = params.RES(i, 1:2);    
    M1_multires{i} = resampleSphericalGrid(squeeze(multiResM1(i, :,:,:)), res(1), 0); 
    M2_multires{i} = resampleSphericalGrid(squeeze(multiResM2(i, :,:,:)), res(1), 0); 
    M3_multires{i} = resampleSphericalGrid(squeeze(multiResM3(i, :,:,:)), res(1), 0); 
      
end

%% Compute the Qmaps at each resolution
Q1_multires = cell(1, params.nres);
Q2_multires = cell(1, params.nres);
Q3_multires = cell(1, params.nres);
d = 3;
for i=1:params.nres,
    res = params.RES(i, 1:2);
    
    Q11 = surface2qnew_qie(mySurf(M1_multires{i}, 0), res, params.mysmall); % , sinTheta, d, N, res);    
    Q1_multires{i} = permute(reshape(Q11, [3, res(1), res(2)]), [2 3 1]);
    
    Q22 = surface2qnew_qie(mySurf(M2_multires{i}, 0), res, params.mysmall); % , sinTheta, d, N, res);    
    Q2_multires{i} = permute(reshape(Q22, [3, res(1), res(2)]), [2 3 1]);
    
    Q33 = surface2qnew_qie(mySurf(M3_multires{i}, 0), res, params.mysmall); % , sinTheta, d, N, res);    
    Q3_multires{i} = permute(reshape(Q33, [3, res(1), res(2)]), [2 3 1]);
    
end

% shooting vector between M1 and M2
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
nsteps = numel(S);
Qpath_multires = cell(1, nsteps);
linearPath = cell(1, nsteps);

for j=1:nsteps
    Q = cell(1, params.nres);
    F = cell(1, params.nres);
    s = S(j);
        
    for i = 1:params.nres,   
        q    = squeeze(Q3_multires{i});
        qV   = QV_multiRes{i};      
        Q{i} = q + s * qV;
        
        f    = squeeze(M3_multires{i});
        fV   = V_multiRes{i};
        F{i} = f + s * fV;
    end
    
    Qpath_multires{j} = Q;
    linearPath{j} = F;
end

%% the inversion of the entire path
geodesicPath{1} = M3_multires;
res = params.RES(end, :);
L   = params.LLS{end};
l   = L(end);
for ii=2:nsteps,
    
    % this is the Q to invert
    Q_multires  = Qpath_multires{ii};
    
      
    %% inversion (comment / uncomment to select the appropriate option)
    % starting with the linear path
    f0_multires = linearPath{ii}; % geodesicPath{ii-1}; %    
    f0  = f0_multires{1};     
    res = params.RES(1, :);
    
    if params.basis_type == 1, % harmonic basis
        load( ['HarmonicBasis/harmonic_basis_Res' num2str(res(1)) '_L' num2str(l) '.mat'], 'B');  
        [cn, f0n] = reconstruct_surfaceR3(mySurf(f0, 0), params.myInner, B, res);
    else % PCA basis
        load( ['PCABasis_HumanShapes/PCABasisHuman_Res' num2str(res(1)) '.mat'], 'B', 'Mu');  
        npca_basis = 100;
        B = B(:,:, 1:npca_basis);
        [cn, f0n] = reconstruct_surfaceR3_PCA(mySurf(f0, 0), B, res, Mu);
        clear('Mu');
    end
    clear('B');
    
 
  
%     Qs   = Q_multires;  %% Qs should be res x res x 3 not 3 x XX
%     LL   = LLS;
%     RESn = RES;
    [f_multires, cn_multires, E, ~] = multiresSRNFInversion(Q_multires,  cn, params); 
    
    
%    invertSRNF(Q_multires, params);
%     multiresSRNFInversion(Qs, ...
%                                                         RESn, LL, ...
%                                                         cn, ...
%                                                         params);
%     % disp(E{end}(end));
    % Use plot(E{end}) to see the evolution of the reconstruction error
    % You will notice that it plateaus quickly - thus you don't 5000
    % iterations
    % 200 are usually sufficient    
    
    % Save
    geodesicPath{ii} = f_multires;
    
     pause(.1);  
end
% geodesicPath{nsteps} = M2_multires;

% the final geodesic at hires
F = [];
for i=1:nsteps,
    R  = geodesicPath{i};
    F(:,:,i,:) = R{end}; % rotateGrid(R, O');
end

% visualize the path
fig_geodesics = 2000;
h1 = figure(fig_geodesics); clf;
DisplayGeod(F, fig_geodesics, 30, 30);

%% Save - uncomment if needed
% save([outDir fname '_geodesic.mat'], 'F', 'M1', 'M2');
% saveas(h1, [outDir fname '_geodesic.fig']);

% the linear path
F = [];
for i=1:nsteps,
    R = linearPath{i}; % Fpath_multires{i};    
    F(:,:,i,:) = R{end}; % rotateGrid(R, O');   
end

% visualize the path
fig_linear = 4000;
h2 = figure(fig_linear); clf;
DisplayGeod(F, fig_linear, 30, 30);

% save([outDir fname '_linear.mat'], 'F', 'M1', 'M2');
% saveas(h2, [outDir fname '_linear.fig']);

