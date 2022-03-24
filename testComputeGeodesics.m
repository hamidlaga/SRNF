%% README
% This script demonstrates how compute geodesics, using SRNF inversion, between two registered surfaces.
% The surfaces need to be pre-registered.
%
%
% The algorithm uses a multiscale representation of surfaces. The examples
% provided here have been generated using sphercial wavelet transform.
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
% The script below takes two surface, computes their SRNFs. It then
% computes a linear path in the SRNF space and inverts it to get the
% geodesic in the original space of surfaces
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
clear all;

configure_paths;

harmonicBasisDir = '../HarmonicBasis/';

dataDir = './sample_surfaces/';
fname   = 'Cat1_1_Cat1_2';
 
%% Load two registered surfaces
load([dataDir fname, '.mat'],  'multiResM1', 'multiResM2');  
M1 = squeeze(multiResM1(end, :,:, :));
M2 = squeeze(multiResM2(end, :,:, :));

%% Parameters
% the resolutions
RES = [128, 128; ... % 25, 25; ...
        8,  8; ...
        16, 16; ...
        25, 25; ...
         50 50;  ...
       ];
   
% frequencies
LLS = { [36]; ...
        [2:8];  ...      
        [8:2:34];    ...   
        [8:2:34]; ... 
        [34,36]; ... 
        };    
nres = min(size(RES, 1), numel(LLS));  

nsteps = 5;           % number of samples along the geodesic
correctSRNFmap   = 0; % set it to 1 to apply the correction suggested by Eric Klassen

srnf_interp_mode = 1; % 0 (default) - use Euclidean distance in SRNF space
                      % set it to 1 if you want to use arc length distance
ix     = 10; % [4, 5, 6, 10]; % the scales that we are going to use (ignored if not starting from a sphere)

% The parameters for SRNF inversion
params.mysmall  = 0.01;
params.d        = 3;
params.myInner  = @innerS2;
params.cutoff   = 1e-6;
params.stepsize = .1;       % 0.01;
params.itermax  =  20000; % 10000;     % Max No. of iterations

params.initMode = 1;        % 0 - Initialize the inversion using a sphere
                            % 1 - Initialize the inversion using a point on
                            %     the linear path
                            % 2 - Initialize the inversion using the
                            %     previously estimated shape
                            % 3 - starts with a given path (see: need to
                            % set params.initPath to the file that contains
                            % the initial geodesic
                            
params.basis_type = 1;      % 1 - Spherical Harmonic basis (default)
                            % 2 - set it to 2 if you want to use PCA basis
                            % (human shapes). In that case you will need to
                            % go to multiresSRNInversion.m and set the path
                            % to the right PCA basis, and eventulally
                            % adjust the paths to Sebastian's/Ian's cods
                            
% params.npca_basis = 20;
% params.PCA_BASIS_FNAME_PREFIX = 'PCABasis_HumanShapes/PCABasisHuman_Res'; % 'PCABasis_cat/PCABasisCat_Res';
                            
toSave            = 1;      % set it to 1 if you want to save the computed gedoesics and figures

%% Normalizing the surfaces
multiResM1(nres, :,:,:) = multiResM1(end, :,:,:);
multiResM2(nres, :,:,:) = multiResM2(end, :,:,:);
multiResM1 = multiResM1(1:nres, :,:,:);
multiResM2 = multiResM2(1:nres, :,:,:);

multiResM1 = normalize_translation_multiresSurf(multiResM1);
multiResM2 = normalize_translation_multiresSurf(multiResM2);

% 
multiResM1 = normalize_scale_multiresSurf(multiResM1);
multiResM2 = normalize_scale_multiresSurf(multiResM2);

%% Perform a rigid registration step (use the qmaps not the f) witht the sphere
[~, O] = rigidAlign(squeeze(multiResM1(end, :,:,:)), squeeze(multiResM2(end, :,:,:)), ...
                    surfaceCodeDir, currDir);

% Apply O to M2
for i=1:size(multiResM2, 1),
    F = squeeze(multiResM2(i, :,:,:));
    multiResM2(i, :,:,:) = rotate3D(F, O);
end

%     0.7638    0.5935    0.2537];
O = eye(3, 3);

% visualize the aligned surfaces
for i=size(multiResM1, 1):size(multiResM1, 1)
    figure(301); clf; hold on;
    F = squeeze(multiResM1(i, :,:,:));
    surface(squeeze(F(:,:, 1)), squeeze(F(:,:, 2)), squeeze(F(:,:, 3)));
    axis equal;
    cameramenu
    %pause
end
for i=size(multiResM1, 1):size(multiResM2, 1)
   % figure(401); clf; hold on;
    F = squeeze(multiResM2(i, :,:,:));
    surface(squeeze(F(:,:, 1)), squeeze(F(:,:, 2)), squeeze(F(:,:, 3)));
    axis equal; cameramenu;
    axis equal;
    cameramenu
   % pause
end
% disp('Press enter to continue ...');
% pause

M1 = squeeze(multiResM1(end, :,:,:));        
M2 = squeeze(multiResM2(end, :,:,:));        

%% Multi-resolution representation of the surface M_orig
for i=1:nres,
    res = RES(i, 1:2);    
    M1_multires{i} = resampleSphericalGrid(squeeze(multiResM1(i, :,:,:)), res(1), 0);  
    M2_multires{i} = resampleSphericalGrid(squeeze(multiResM2(i, :,:,:)), res(1), 0);        
end

% maybe not needed
M1_multires = normalize_scale_multiresSurf(M1_multires);
M2_multires = normalize_scale_multiresSurf(M2_multires);

%% Compute the Qmaps at each resolution
Q1_multires = cell(1, nres);
Q2_multires = cell(1, nres);

for i=1:nres,

    % Using Qian's code
    res = RES(i, :);
    mySmall = 0.01;
    Q11 = surface2qnew_qie(mySurf(M1_multires{i}, 0), res, mySmall); % , sinTheta, d, N, res);    
    Q1_multires{i} = permute(reshape(Q11, [3, res(1), res(2)]), [2 3 1]);
    
    Q22 = surface2qnew_qie(mySurf(M2_multires{i}, 0), res, mySmall); % , sinTheta, d, N, res);
    Q2_multires{i} = permute(reshape(Q22, [3, res(1), res(2)]), [2 3 1]);
    
end


%% Compute the linear path between Q1 and Q2
% linear path
linearPath = cell(1, nsteps);

for j=1:nsteps
    F = cell(1, nres);
    s = (j-1)/(nsteps-1);
    for i = 1:nres,   
        f1   = M1_multires{i};
        f2   = M2_multires{i};
        F{i} = f1*(1-s)+ f2*s;
    end
    
    linearPath{j} = F;
end

%%
Qpath_multires = cell(1, nsteps);

for j=1:nsteps
    Q = cell(1, nres);
    s = (j-1)/(nsteps-1);
       
    for i = 1:nres,   

        q1  = Q1_multires{i};
        q2  = Q2_multires{i};
        
        if srnf_interp_mode == 1,
            % nicer interpolation of SRNFs
            Q{i} = transportSRNFdeformation(q1, q2, s, q1);
        else 
            % linear interpolation of SRNFs (default)
            Q{i} = q1*(1-s) + q2*s;
        end

        % Correcting the SRNF map (Eric's idea)
        if correctSRNFmap,
            [Q{i},  intQ,  norm_intQ, intQn, norm_intQn] = correctQmap(Q{i}); 
        end
    end
    
    Qpath_multires{j} = Q;
end

%% the inversion of the entire path
geodesicPath{1} = M1_multires;
% res = RES(end, :);  
% L   = LLS{size(RES, 1)};
% l   = L(end);
% load( ['HarmonicBasis/harmonic_basis_Res' num2str(res(1)) '_L' num2str(l) '.mat']);     

cn_multires = [];
for ii = 2:nsteps-1,
    
    % this is the Q to invert
    Q_multires =  Qpath_multires{ii};
    
    
    if params.initMode  == 0,   % initialize with a sphere
        
        [f_multires, cn_multires, E] = multiresSRNInversion(Q_multires, ...
                                                        RES(1:nres, :), LLS, ...
                                                        [], ...
                                                        params);
        
    elseif params.initMode == 1,
        
        res = RES(nres, :);
        L   = LLS{nres};
        l   = L(end);
        
        if params.initMode  == 1,       % linear path
            f0_multires =   linearPath{ii};
            
        elseif params.initMode  == 2,   % previously estimated shape
            f0_multires = geodesicPath{ii-1}; % 
        end
        
        f0  = f0_multires{end};         % the initial surface
        
        % project it to the harmonic basis
        if params.basis_type == 1,  % spharmonics
            load( [harmonBasisDir 'harmonic_basis_Res' num2str(res(1)) '_L' num2str(l) '.mat'], 'B');  
            [cn, ~] = reconstruct_surfaceR3(mySurf(f0, 0), params.myInner, B, res);
            clear('B');
        else
            
            % load PCA basis
%             load( [params.PCA_BASIS_FNAME_PREFIX num2str(res(1)) '.mat']);
%             [cn, ~] = reconstruct_surfaceR3_PCA(mySurf(f0, 0), B, res, Mu);
            cn = []; 
        end
        
        % inversion
        Qs{1} = Q_multires{nres};
        LL{1} = LLS{nres};
        [f_multires, cn_multires, E] = multiresSRNInversion(Qs, ...
                                                            RES(1:nres, :), LL, ...
                                                            cn, ...
                                                            params);
    else % take a pre-computed path
        
        res = RES(nres, :);
        L   = LLS{nres};
        l   = L(end);
        
        load(params.initPath, 'F');
        f0 = resampleSphericalGrid( squeeze(F(:,:, 4, :)), res(1), 0 ); 
        
        q0 = surface2qnew_qie(mySurf(f0, 0), res, .01); %  mySmall); % , sinTheta, d, N, res);
        q0 = permute(reshape(q0, [3, res(1), res(2)]), [2 3 1]);
        
        load( [harmonBasisDir 'harmonic_basis_Res' num2str(res(1)) '_L' num2str(l) '.mat'], 'B');  
        [cn, ~] = reconstruct_surfaceR3(mySurf(f0, 0), params.myInner, B, res);
        clear('B');
        
        % inversion
        Qs{1} = Q_multires{nres};
        
        Qs{1} = correctOrientationOfSRNF(Qs{1}, q0);
        
        visualizeNormalField(f0, Qs{1}); %q1, q2,  Qs{1});
        
        LL{1} = LLS{nres};
        [f_multires, cn_multires, E] = multiresSRNInversion(Qs, ...
                                                            RES(1:nres, :), LL, ...
                                                            cn, ...
                                                            params);
        
    end
    
    % the final energy
    disp(E{end}(end));
    
    % save
    geodesicPath{ii} = f_multires;
    
    %% visualization
%     for i=nres:nres,
%         f4 = f_multires{1}; %{i};
%         res= RES(i, :);
% 
%         %% visualize the reconstructed surface
% %         M4 = reshape(f4,[d,res]);
% %         M4 = permute(M4,[2,3,1]);
% 
%         h1 = figure(100); clf; 
%         dispSurfR3(M1_multires{i}, res, 3);  cameramenu
% 
%         h2 = figure(200); clf;
%         dispSurfR3(f4,res,3);cameramenu
% 
%     end
   % disp('Press anykey to continue ...');
    pause(.1);
    
end
geodesicPath{nsteps} = M2_multires;

%% Recover the gedoesic at high res (ignore for the moment - it's not working)
% [res_n(1), res_n(2), ~] = size(squeeze(multiResM1(end, :,:,:))); %  = [200, 200];
% % res_n = [50, 50]; 
% M1 = resampleSphericalGrid(squeeze(multiResM1(end, :,:,:)), res_n(1), 0); 
% M2 = resampleSphericalGrid(squeeze(multiResM2(end, :,:,:)), res_n(1), 0);
% Fgeod_hires = [];
% Fgeod_hires(:,:, 1, :) = M1;
% Fgeod_hires(:,:, nsteps, :) = M2;
% 
% Fref = M1;
% 
% for i=1:nsteps,
%     R = geodesicPath{i};
%     R = R{end};
%     f0 = reshape(R, [d,res]);
%     F  = permute(f0,[2,3,1]);       % the intermediate shape
%     
%     Fgeod_hires(:,:, i, :) = upscale(F, Fref);
% end
% fig_geodesics_hires = 5000;
% hgeod_hires = figure(fig_geodesics_hires); clf;
% DisplayGeod(Fgeod_hires, fig_geodesics_hires, 30, 30);


%% Visualizing the computed geodesics
% the final geodesic at hires
Fgeod = [];
for i=1: nsteps,
    R  = normalize_scale_multiresSurf(geodesicPath{i});    
    Fgeod(:,:,i,:) = rotateGrid(R{end}, inv(O));
end

fig_geodesics = 2000;
hgeod = figure(fig_geodesics); clf;
DisplayGeod(Fgeod, fig_geodesics,30,30);

% the linear path at hires
M1 =  M1_multires{nres}; 
M2 =  M2_multires{nres}; 
Flin = [];
for i=1:nsteps,
    s = (i-1)/(nsteps-1);
    Flin(:,:,i,:) = rotateGrid(M1*(1-s)+ M2*s, inv(O));    
end

fig_linear = 4000;
hlin = figure(fig_linear); clf;
DisplayGeod(Flin, fig_linear,30,30);

%% Saving results and figures
if toSave,
    F = Fgeod;
    save([outDir fname_out '_geodesic.mat'], 'F');
    saveas(hgeod, [outDir fname_out '_geodesic.fig']);

    F = Flin;
    save([outDir fname_out '_linearpath.mat'], 'F');
    saveas(hlin, [outDir fname_out '_linearpath.fig']);
end

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



