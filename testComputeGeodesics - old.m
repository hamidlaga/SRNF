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

%% Multi-resolution representation of the surface M_orig
% the resolutions
RES  = [8, 8; ...
        16, 16; ...
        25, 25; ...
        50, 50]; 
nres = size(RES, 1);

% The scales (frequencies)
LLS{1} = [2:8]; 
LLS{2} = [8:12];
LLS{3} = [16:2:22];
LLS{4} = [22]; % [22];
% LLS{1} = [22]; % [22];

multiResM1(nres, :,:,:) = multiResM1(end, :,:,:);
multiResM2(nres, :,:,:) = multiResM2(end, :,:,:);
M1_multires = cell(1, nres);
M2_multires = cell(1, nres);

for i=1:nres,
    res = RES(i, 1:2);    
    M1_multires{i} = mySurf(resampleSphericalGrid(squeeze(multiResM1(i, :,:,:)), res(1), 0), 0); 
    M2_multires{i} = mySurf(resampleSphericalGrid(squeeze(multiResM2(i, :,:,:)), res(1), 0), 0); 

    % visualization just for test (uncomment if needed)
%     figure(i); clf;
%     dispSurfR3(M1_multires{i}, res, 3);
%     cameramenu;
%     
%     figure(i*1000); clf;
%     dispSurfR3(M2_multires{i}, res, 3);
%     cameramenu;
% 
%     pause

end

%% Compute the Qmaps at each resolution
Q1_multires = cell(1, nres);
Q2_multires = cell(1, nres);
for i=1:nres,
    res = RES(i, 1:2);
    
    [Theta,~] = genGridSphr(res);
    mysmall = 0;
    sinPhi  = genSinPhi(Theta, mysmall);
    N       = prod(res);
    d       = 3;
    
    [Mu, Mv] = partialF(M1_multires{i}, sinPhi, d, N, res);
    [Q1_multires{i},  ~, ~, ~] = map2qnew(Mu, Mv, N);

    [Mu, Mv] = partialF(M2_multires{i}, sinPhi, d, N, res);
    [Q2_multires{i},  ~, ~, ~] = map2qnew(Mu, Mv, N);

end

%% Compute the linear path between Q1 and Q2
nsteps = 3;
Qpath_multires = cell(1, nsteps);

for j=1:nsteps
    Q = cell(1, nres);
    s = (j-1)/(nsteps-1);
    for i = 1:nres,   

        q1  = squeeze(Q1_multires{i});
        q2  = squeeze(Q2_multires{i});
        
        Q{i} = q1*(1-s)+ q2*s;
    end
    
    Qpath_multires{j} = Q;

end

% linear path
linearPath = cell(1, nsteps);

for j=1:nsteps
    F = cell(1, nres);
    s = (j-1)/(nsteps-1);
    for i = 1:nres,   

        f1  = squeeze(M1_multires{i});
        f2  = squeeze(M2_multires{i});
        
        F{i} = f1*(1-s)+ f2*s;
    end
    
    linearPath{j} = F;
end

%% the inversion of the entire path
params.mysmall  = 0;
params.d        = 3;
params.myInner  = @innerS2;
params.cutoff   = 1e-6;
params.stepsize = .1;       %0.01;
params.itermax  = 10000; 

geodesicPath{1} = M1_multires;
res = RES(end, :);
L   = LLS{end};
l   = L(end);
 
cn_multires = [];
for ii=2:nsteps-1,
    
    % this is the Q to invert
    Q_multires =  Qpath_multires{ii};
    
    % use the previous surface to initialize the inversion
    f0_multires  =  linearPath{ii}; 
    f0  = f0_multires{end};
    res = RES(end, :);

    if isempty(cn_multires)
        load( [ harmonicBasisDir 'harmonic_basis_Res' num2str(res(1)) '_L' num2str(l) '.mat'], 'B');  
        [cn, ~] = reconstruct_surfaceR3(f0, params.myInner, B, res);
        clear('B');
    else
        cn = cn_multires{end};
    end
    
    % inversion
    [f_multires, cn_multires, E] = multiresSRNInversion(Q_multires, RES, LLS, [], params);
    disp(E(end));

    % save
    geodesicPath{ii} = f_multires;
    
    %% visualization
    for i=1:nres,
        f4 = f_multires{i};
        res= RES(i, :);

        %% visualize the reconstructed surface
%         M4 = reshape(f4,[d,res]);
%         M4 = permute(M4,[2,3,1]);

        h1 = figure(100); clf; 
        dispSurfR3(M1_multires{i}, res, 3);  cameramenu

        h2 = figure(200); clf;
        dispSurfR3(f4,res,3);cameramenu

    end
    disp('Press anykey to continue ...');
    pause (.1);
    
end
geodesicPath{nsteps} = M2_multires;

% the final geodesic at hires
F = [];
for i=1:nsteps,
    R = geodesicPath{i};
    R = R{end};
    f0 = reshape(R, [d,res]);
    R  = permute(f0,[2,3,1]);

    F(:,:,i,:) = R;
end

% visualize the path
fig_geodesics = 2000;
h = figure(fig_geodesics); clf;
DisplayGeod(F, fig_geodesics,30,30);

% save([outDir fname '_geodesic.mat'], 'F');
% saveas(h, [outDir fname '_geodesic.fig']);

% the linear path
M1 = M1_multires{end};
f0 = reshape(M1, [3,res]);
M1  = permute(f0,[2,3,1]);

M2 = M2_multires{end};
f0 = reshape(M2, [3,res]);
M2  = permute(f0,[2,3,1]);

F = [];
for i=1:nsteps,
    s = (i-1)/(nsteps-1);
    F(:,:,i,:) = M1*(1-s)+ M2*s;
    
end

% visualize the path
fig_geodesics = 4000;
h = figure(fig_geodesics); clf;
DisplayGeod(F, fig_geodesics,30,30);



