function f_reconstructed_multires = invertSRNF_PCABasis(Q_multires, params)

if nargin < 2,
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
    
end

%% Now start the inversion process (the input is Q_multires and RES, the resolutions)
cn = [];    % initial surface
tic
[f_reconstructed_multires, cn, E, dE] = multiresSRNFInversion(Q_multires, cn, params);
toc