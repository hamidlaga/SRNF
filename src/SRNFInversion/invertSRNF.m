function f_reconstructed_multires = invertSRNF(Q_multires, params)

if nargin < 2,
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

    params.npca_basis = 20;
    params.PCA_BASIS_FNAME_PREFIX = 'PCABasis_HumanShapes/PCABasisHuman_Res'; % 'PCABasis_cat/PCABasisCat_Res';
    params.harmonic_basis_homedir = '/Users/hamidlaga/Dropbox/Home/Applications/3DShapeStatistics/Release/InverseQmap/HarmonicBasis/';                            
    params.toSave = 1;      % set it to 1 if you want to save the computed gedoesics and figures    
    
    params.cn = [];         % initial surface - sphere
    
end

%% Now start the inversion process (the input is Q_multires and params.RES, the resolutions)
cn = params.cn;    % initial surface
tic
[f_reconstructed_multires, cn, E, dE] = multiresSRNFInversion(Q_multires,  cn, params); 
toc