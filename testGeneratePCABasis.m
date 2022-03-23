%% Make PCA basis

% Computes linear statistics from the registered 3D shapes
%

clear all;
close all;
imtool close all;

drive = 'F:/';      %PC

homeDir        = [drive 'Home/Hamid/Work/Applications/3DShapeStatistics/code/'];        
surfaceCodeDir = [homeDir 'Surface/'];
currDir        = pwd;

dataDir = [homeDir '../3DModels/MPI-EG2009/multires_surfaces/']; % [homeDir '../3DModels/TOSCA/multi_res/'];
outDir  = './PCABasis_HumanShapesFAUST/';  % './PCABasis_cat/';  %PCABasisHuman_Res';
mkdir (outDir);
outfname =  [outDir 'PCABasisHuman_Res'];

D       = dir([dataDir '*.mat']); %  'cat*.mat']);

refIx   = 1;
res     = [32, 32]; % [32, 32]; %  [64, 64]; % [32, 32]; % desired resolution

nObjects = numel(D);

%% Load the shapes
load([dataDir D(1).name], 'multiResM');

M = resampleSphericalGrid(squeeze(multiResM(end, :,:,:)), res(1), 0);
[Theta, Phi]      = genGridSphr(res);
[A1, multfact1, ~] = area_surf_closed(M);
M        = center_surface_closed(M, multfact1,A1,Theta);
[A1, ~, ~] = area_surf_closed(M);
M = M /sqrt(A1);

% muF = M;
P1  = M;

[N,~,~] = size(M);

P1 = reshape(P1, N*N, []);
P1 = P1';

Shapes = [];

figure; clf; hold on;
for i=1:nObjects,
 
    disp([num2str(i) '/' num2str(nObjects) ' ...']);
    load([dataDir D(i).name]);
    M = resampleSphericalGrid(squeeze(multiResM(end, :,:,:)), res(1), 0);
    [Theta,Phi] = genGridSphr(res);
    [A1,multfact1,~] = area_surf_closed(M);
    M = center_surface_closed(M, multfact1,A1,Theta);
    [A1,~,~] = area_surf_closed(M);
    M = M /sqrt(A1);
    
    %rigidly align M to S0
    Q1 = reshape(M, N*N, []); 
    Q1 = Q1';
    A  = P1 * Q1';
    [U,S,V] = svd(A); 
    O       = U*V';         % aligns Q on P
    S1      = rotate3D(M, O);
    
    Shapes = [Shapes S1(:)];    
end

%% Computing the linear statistics
nObjects = size(Shapes, 2);
nModes   = nObjects-1;
[Mu, eigenVectors, EVals] = performEigenAnalysis(Shapes(:, 1:end), nModes); % , muF(:));
muF = reshape(Mu, [N, N, 3]);

%% visualization
cmap = hsv;
cmap = [cmap; 1 0.95 0.9; 0.7 0.7 0.7];
ncolors = size(cmap, 1);
[h, w, ~] = size(muF);

% the mean without landmarks
set(0,'DefaultFigurePosition',[200 300 560 420]);
hfig_mean = 10;
figure(hfig_mean), clf;
surface(squeeze(muF(:, :, 1)), squeeze(muF(:, :, 2)), squeeze(muF(:, :, 3)), ...
         ncolors * ones( h, w), 'SpecularExponent',100);
axis off equal;
cameramenu;

lighting phong
camlight('headlight')
shading 'flat';
colormap(cmap);
set(gcf, 'Renderer', 'zbuffer');
set(gcf, 'Color', [1, 1, 1]);


%% doing some reconstruction
% S = Shapes(:, 89);
% Sn = S - Mu;
% nModes = size(eigenVectors, 2);
% 
% % projection
% Cn = zeros(1, nModes);
% for i=1:nModes
%    
%    V = eigenVectors(:, i);
%    Cn(i) = dot(Sn,  V);
% end
% 
% Sn = Mu + eigenVectors(:, 1:nModes) * Cn(:);
% 
% % visualization
% M1 = reshape(S, [N, N, 3]);
% M2 = reshape(Sn, [N, N, 3]);
% 
% figure(10), clf;
% surface(squeeze(M1(:, :, 1)), squeeze(M1(:, :, 2)), squeeze(M1(:, :, 3)), ...
%          ncolors * ones( N, N), 'SpecularExponent',100);
% axis off equal;
% cameramenu;
% 
% figure(2), clf;
% surface(squeeze(M2(:, :, 1)), squeeze(M2(:, :, 2)), squeeze(M2(:, :, 3))+.1, ...
%          ncolors * ones( N, N), 'SpecularExponent',100);
% axis off equal;
% cameramenu;
% pause

% the modes of variation
for modeId = 1:min(size(eigenVectors, 2), 10),  % Let's visualize 10 modes

    % close all;
    eVal  = sqrt(EVals(modeId));        % standard deviation

    f = -2:0.5:2; % [-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3];
    Variations =  f* eVal;

    h = figure(modeId * 1000);
    clf;
    hold on;
    colors = {'.r' ; '.g' ; '.b'; 'xg'; 'xr'};
    
    shift_x = 1;
    for i=1:length(Variations)
        sc = Variations(i);
        S  = reshape(Mu + sc * eigenVectors(:, modeId), [N, N, 3] ); 

        surface(squeeze(S(:, :, 1)) + shift_x * i, squeeze(S(:, :, 2)), squeeze(S(:, :, 3)));
    end   
  
    axis equal;
    cameramenu;
    axis off;
    
    lighting phong;
    light('Position',[0 0 1],'Style','local')
%    camlight('headlight');
    shading 'flat';
    set(gcf, 'Renderer', 'zbuffer');
    
	% savefig(h, [outDir 'Hamate_mode_' num2str(modeId) '.fig']);
    disp('Press any key to visualize the next mode of variation ...');
    pause; % (.1)
    
end

%% saving the PCA basis
% B = permute(eigenVectors, [2, 1]);
B = reshape(eigenVectors, [N, N, 3, size(eigenVectors, 2)]);
B = permute(B, [ 3 1 2 4]);
B = reshape(B, [3, N*N, size(B, 4)]);
Mu = permute(muF, [3, 1 , 2]);
Mu = reshape(Mu, 3, []);

%% generating the derivatives of the basis
Bn = zeros(size(B, 1), size(B, 2), size(B, 3)+1);
Bn(:,:, 1) = Mu;
Bn(:,:, 2:end) = B;

[Bnu,Bnv] = BasisDiffR3(Bn,[N N]);
Muu = squeeze(Bnu(:,:, 1));
Muv = squeeze(Bnv(:,:, 1));
clear Bn;

Bu = Bnu(:,:, 2:end); clear Bnu;
Bv = Bnv(:,:, 2:end);


save([outfname num2str(N) '.mat'], ...
    'B', 'Bu', 'Bv', 'Mu', 'Muu', 'Muv');




