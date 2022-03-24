%% compute geodesics by SRNF inversion of the linear path

%
% Uses multires surfaces generated with sph wavelet decomposition
%
%
clear all;
close all;

configure_paths;

mysmall = 0.01;

dataDir = './sample_surfaces/';
outDir  = './output/'; 
homeDir        = [pwd '/']; 
surfaceCodeDir = [homeDir '../Surface/'];  %% Specify the path to Surface code
currDir        = pwd;

%% The shapes
fname1 = 'dfaust_0';
fname2 = 'dfaust_1'; 

load([dataDir fname1, '.mat'], 'multiResM');  multiResM1 = multiResM;
load([dataDir fname2, '.mat'], 'multiResM');  multiResM2 = multiResM;

% normalize for translation
multiResM1 = normalize_translation_multiresSurf(multiResM1);
multiResM2 = normalize_translation_multiresSurf(multiResM2);

% normalize for scale
multiResM1 = normalize_scale_multiresSurf(multiResM1);
multiResM2 = normalize_scale_multiresSurf(multiResM2);

%% remove the poles from narrow features (this is done manually at this stage)
O = eye(3, 3);

M1 = squeeze(multiResM1(end, :,:,:));        
M2 = squeeze(multiResM2(end, :,:,:));        

%% test here the optimal rotation of the grid for M2
% the standard sphere
[res(1), res(2), ~] = size(M1);
[M2n, O] = alignSurfaceGridWithSphere(M2, surfaceCodeDir, currDir); 
F = M2n;

figure(10); clf;
surface(F(:,:,1),F(:,:,2),F(:,:,3));
axis equal;
cameramenu;
pause

% O = getRandomRotationMatrix;
for i=1:size(multiResM1, 1),
    M1 = squeeze(multiResM1(i, :,:,:));     
    multiResM1(i, :,:,:) = rotateGrid(M1, O);
    
    M2 = squeeze(multiResM2(i, :,:,:));     
    multiResM2(i, :,:,:) = rotateGrid(M2, O);
end

%% the resolutions
RES = [64, 64]; 

%% The scales (frequencies)
LLS{1} = [36]; 

nres = min([numel(LLS), size(RES, 1)]);

multiResM1(nres, :,:,:) = multiResM1(end, :,:,:);
multiResM2(nres, :,:,:) = multiResM2(end, :,:,:);

M1_multires = cell(1, nres);
M2_multires = cell(1, nres);

for i=1:nres,
    res = RES(i, 1:2);  
    
    K = resampleSphericalGrid(squeeze(multiResM1(i, :,:,:)), res(1), 0);
    [A1,~,~] = area_surf_closed(K);   A1 = 1; 
    M1_multires{i} = mySurf(K / sqrt(A1), 0); 
    
    K = resampleSphericalGrid(squeeze(multiResM2(i, :,:,:)), res(1), 0);
    [A2,~,~] = area_surf_closed(K); A1 = 1;
    M2_multires{i} = mySurf(K / sqrt(A2), 0);     
     
%    % visualization just for test
    figure(i); clf;
    dispSurfR3(M1_multires{i}, res, 3);
    cameramenu;
    hold on;
    % figure(i*1000); clf;
    dispSurfR3(M2_multires{i}, res, 3);
    cameramenu;

    % pause

end

%% Compute the Qmaps at each resolution
Q1_multires = cell(1, nres);
Q2_multires = cell(1, nres);
for i=1:nres,
    res = RES(i, 1:2);    
    [Theta,Phi] = genGridSphr(res, mysmall);    
    sinTheta  = genSinPhi(Theta, mysmall);
    N       = prod(res);
    d       = 3;
    
    [Mu, Mv] = partialF(M1_multires{i}, sinTheta, d, N, res);
    [Q1_multires{i},  ~, ~, ~] = map2qnew(Mu, Mv, N);    
    
    [Mu, Mv] = partialF(M2_multires{i}, sinTheta, d, N, res);
    [Q2_multires{i},  ~, ~, ~] = map2qnew(Mu, Mv, N);
    
end

%% Compute the linear path between Q1 and Q2
nsteps = 3; % just the mean

% linear path
linearPath = cell(1, nsteps);
for j=1:nsteps
    F = cell(1, nres);
    s = (j-1)/(nsteps-1);
    for i = 1:nres,   
        f1   = squeeze(M1_multires{i});
        f2   = squeeze(M2_multires{i});
        F{i} = f1*(1-s)+ f2*s;
    end
    
    linearPath{j} = F;
end

%% Linear path in the space of SRNfs
Qpath_multires = cell(1, nsteps);
Rots_multires  = cell(1, nsteps); % holds the rotations that have been applied to the grid
                                  % apply their inverse to get the original
                                  % parameterizations
for j=1:nsteps
    Q = cell(1, nres);
    Rots = cell(1, nres);
    s = (j-1)/(nsteps-1);
    
    for i = 1:nres,   
        res = RES(i, 1:2);    
        [Theta,Phi] = genGridSphr(res, mysmall);    
        sinTheta    = genSinPhi(Theta, mysmall);
        dphi    = Phi(2, 1) - Phi(1, 1);     % 2*pi/(n);
        dtheta  = Theta(1, 2) - Theta(1, 1);
        N       = prod(res);
    
        q1  = squeeze(Q1_multires{i});
        q2  = squeeze(Q2_multires{i});
        Q{i} = q1*(1-s) + q2*s;  %    interpolateSRNF(q1, q2, s); %       
        
    end
    
    Qpath_multires{j} = Q;
end

%% the inversion of the entire path
params.mysmall  = 0;
params.d        = 3;
params.myInner  = @innerS2;
params.cutoff   = 1e-6;
params.stepsize = .1; %.5;       %0.01;
params.itermax  = 15000; 
params.basis_type = 2; % 2 for PCA, 1 for harmonics

geodesicPath{1}      = M1_multires;
geodesicPath{nsteps} = M2_multires;
res = RES(end, :);
L   = LLS{end};
l   = L(end);
% load( ['HarmonicBasis/harmonic_basis_Res' num2str(res(1)) '_L' num2str(l) '.mat']);     

cn_multires = [];
Elinpath = zeros(nsteps,1);         % energy of the linear path
Egeodpath = zeros(nsteps,1);        % energy of the geodesic path
for ii = 2:nsteps-1,
    
    % this is the Q to invert
    Q_multires  = Qpath_multires{ii};    
      
    %% inversion (comment / uncomment to select the appropriate option)
    % starting with a sphere 
%     [f_multires, cn_multires, E] = multiresSRNInversion(Q_multires, ...
%                                                         RES, LLS, ...
%                                                         [], ...
%                                                         params);
        
    % starting with the linear path
    f0_multires = linearPath{ii}; % geodesicPath{ii-1}; %    
    f0  = f0_multires{1};   
    figure, dispSurfR3(f0, res, 3);
    pause
    res = RES(1, :);
    
    if params.basis_type == 1, % harmonic basis
        load( [homeDir '/HarmonicBasis/harmonic_basis_Res' num2str(res(1)) '_L' num2str(l) '.mat'], 'B');  
        [cn, f0n] = reconstruct_surfaceR3(f0, params.myInner, B, res);
    else % PCA basis
        load( [homeDir 'PCABasis/DFAUST_PCA_with_derivatives_' num2str(res(1)) '.mat'], 'B', 'Mu');  
        [cn, f0n] = reconstruct_surfaceR3_PCA(f0, B, res, Mu);
        clear('Mu');
    end
    clear('B');
    
    Qs = Q_multires;
    LL = LLS;
    RESn = RES;
    [f_multires, cn_multires, E] = multiresSRNInversion(Qs, ...
                                                        RESn, LL, ...
                                                        cn, ...
                                                        params);
    
    Elinpath(ii)  = E{1}(1);
    Egeodpath(ii) = E{end}(end);
    
    disp([Elinpath(ii), Egeodpath(ii)]);
    % Save
    geodesicPath{ii} = f_multires;
    
    %% visualization
    for i=nres:nres,
        f4 = f_multires{i};%{i}; % {i}; % {i};
        res= RES(i, :);

        %% visualize the reconstructed surface
        f0 = f0_multires{i};  
        h1 = figure(100); clf; hold on;
        dispSurfR3(f0, res, 3); % M1_multires{nres}
        axis on; hold on;

        % h2 = figure(200); clf;
        dispSurfR3(f4,res,3);
        axis on;

    end
    disp('Press anykey to continue ...');
    pause(.1);
    
end

% the final geodesic at hires
F = [];
for i=1:nsteps,
    R = geodesicPath{i};
    R = R{end};
    f0 = reshape(R, [d,res]);
    R  = permute(f0,[2,3,1]);

    [A1,~,~] = area_surf_closed(R);    
    R =  R/ sqrt(A1); 
    [A1,~,~] = area_surf_closed(R);    
    disp(A1)
    
    F(:,:,i,:) = R;
end

% visualize the path
fig_geodesics = 3000;
h = figure(fig_geodesics); clf;
DisplayGeod(F, fig_geodesics,30,30);
camva(6.0417);
campos([-62.9447    6.7957   -2.27554]);
camup([ 0.1046    0.9944    0.0118]);
% 
% save([outDir fname '_geodesic_srnf_inversion.mat'], 'F', 'Egeodpath');
% saveas(h, [outDir fname '_geodesic_srnf_inversion.fig']);
% print([outDir fname '_geodesic_srnf_inversion.eps'], '-depsc');

% the linear path
M1 = M1_multires{nres};
M2 = M2_multires{nres};
% load( ['HarmonicBasis/harmonic_basis_Res' num2str(res(1)) '_L' num2str(LLS{end}(end)) '.mat']); 
% [~, M1] = reconstruct_surfaceR3( M1, @innerS2, B, res);
% [~, M2] = reconstruct_surfaceR3( M2, @innerS2, B, res);

M1 =  reshape(M1, [3,res]);
M1  = permute(M1,[2,3,1]);

M2 =  reshape(M2, [3,res]);
M2  = permute(M2, [2,3,1]);

F = [];
for i=1:nsteps,
    s = (i-1)/(nsteps-1);
    F(:,:,i,:) = M1*(1-s)+ M2*s;    
end

% visualize the path
fig_geodesics = 4000;
h = figure(fig_geodesics); clf;
DisplayGeod(F, fig_geodesics,30,30);
camva(6.0417);
campos([-62.9447    6.7957   -2.2755]);
camup([ 0.1046    0.9944    0.0118]);

% 
% save([outDir fname '_linearpath.mat'], 'F', 'Elinpath');
% saveas(h, [outDir fname '_linearpath.fig']);
% print([outDir fname '_linpath.eps'],'-depsc');

%% plot the energy of the two paths
% fig_energy = 5000;
% h = figure(fig_energy); clf;
% hold on;
% plot(Egeodpath, '-xr', 'LineWidth', 2);
% plot(Elinpath, '-ob', 'LineWidth', 2);
% saveas(h, [outDir fname '_pathenergies.fig']);
% print([outDir fname '_pathenergies.eps'],'-depsc');

% 
% % the energy
% h =figure;
% plot(E);
% saveas(h, [outDir fname '_energy.fig']);

%% visualization
% for i=1:nres,
%     f4 = f4_multires{i};
%     res= RES(i, :);
%     
%     %% visualize the reconstructed surface
%     M4 = reshape(f4,[d,res]);
%     M4 = permute(M4,[2,3,1]);
%     
%     h1 = figure(100); clf; 
%     dispSurfR3(M_multires{i}, res, 3);  cameramenu
%    
%     h2 = figure(200); clf;
%     dispSurfR3(f4,res,3);cameramenu
%     
%     dirname = ['./report/spherical_wavelet_results/' fname '_multires'];
%     mkdir (dirname);
%     saveas(h1, [dirname '/' fname '_' num2str(res(1)) 'x' num2str(res(2)) '_original.fig']);
%     saveas(h2, [dirname '/' fname '_' num2str(res(1)) 'x' num2str(res(2)) 'reconstructed.fig']);
%    
%     
%     disp('Press anykey to continue ...');
%     pause; % (.1);
% 
% end


