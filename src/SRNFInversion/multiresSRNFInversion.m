function [f4_multires, cn_multires, E, dE] = multiresSRNFInversion(Q_multires0, cn, params)

if ~isfield(params, 'LLS')
    for i=1:size(params.RES, 1)
       params.LLS{i} = [36]; 
    end
end

surfaceCodeDir = params.surfaceCodeDir;

nres = numel(Q_multires0);
Q_multires = cell(nres, 1);

for i=1:nres,
    Q_multires{i} = mySurf(Q_multires0{i}, 0);
end


mysmall = params.mysmall;
d       = params.d;
myInner = params.myInner ;
cutoff  = params.cutoff  ;
stepsize= params.stepsize;
itermax = params.itermax ; 

if ~isfield(params, 'basis_type'),
    params.basis_type = 1;  % sph harmonics by default
end

nres    = min(numel(Q_multires),  size(params.RES, 1));
cn_multires = [];

if params.basis_type == 2,  % PCA basis
    cn = [1; cn]; 
end


for i=1:nres,
    
    disp(['Processing resolution ' num2str(i) ' ...']);    
    res = params.RES(i, 1:2);
    
    [Theta, ~] = genGridSphr(res);    
    sinTheta   = genSinPhi(Theta, mysmall);
    
    q0 = Q_multires{i};         % the q to invert
    
    for l = params.LLS{i}, 
        % load the harmonic basis
        if params.basis_type == 1,
            load( [params.harmonic_basis_homedir 'harmonic_basis_Res' num2str(res(1)) '_L' num2str(l) '.mat']);     
            O  = eye(3, 3);            
            % reshape q0
            q2 = reshape(q0, [3, res]);
            q2 = permute(q0, [2,3,1]);
        else
            % PCA basis
            load( [params.PCA_BASIS_FNAME_PREFIX num2str(res(1)) '.mat']); 
            
            npca_basis =  min(params.npca_basis, size(B, 3));

            B  = B (:,:, 1:npca_basis);
            Bu = Bu(:,:, 1:npca_basis);
            Bv = Bv(:,:, 1:npca_basis);
            
            Bn = zeros(3, size(B, 2), size(B, 3)+1);
            Bn(:,:, 1)     = Mu;
            Bn(:,:, 2:end) = B;
            B = Bn;            
            
            Bnu = zeros(3, size(Bu, 2), size(Bu, 3)+1);
            Bnu(:,:, 1)     = Muu;
            Bnu(:,:, 2:end) = Bu;
            Bu              = Bnu;
            
            Bnv = zeros(3, size(Bv, 2), size(Bv, 3)+1);
            Bnv(:,:, 1)     = Muv;
            Bnv(:,:, 2:end) = Bv;
            Bv = Bnv;            
            
            clear('Muu');
            clear('Muv');          
            
            %% align q0 with the mean shape
            sinTheta  = genSinPhi(Theta, mysmall);
            N       = prod(res);
            d       = 3;

            [mu, mv] = partialF(Mu, sinTheta, 3, prod(res), res);
            [q1, ~, ~, ~] = map2qnew(mu, mv, prod(res));
            q1 = reshape(q1, [3, res]);
            q1 = permute(q1,[2,3,1]);
    
            currDir = pwd;
            
            cd([surfaceCodeDir   'Matlab/ClosedIan']); % 'Mex\ClosedIan']);  % 'Matlab\ClosedIan']); % 
            addpath([surfaceCodeDir 'Matlab/ClosedIan']);
            addpath([surfaceCodeDir 'Matlab/Closed']);

            q2 = reshape(q0, [3, res]);
            q2 = permute(q2,[2,3,1]);
            
            % Rotation that aligns q2 to q1
            O  = optimal_rot_surf_closed(q1, q2, Theta);    
            O  = real(O);
            cd(currDir);
        end
            
        % Align q2 to q1 (either sphere or mean shape)
        q0 = rotate3D(q0, O);
        % 1. inverse mapping, codes in cpp new
        if ~isempty(cn),
            % padding cn with zeros
            if numel(cn) < size(B, 3), % pad with zeros
                cn = [cn; zeros(size(B, 3) - length(cn), 1)];
            else % truncate
                cn = cn(1:size(B, 3)); 
            end
        end
                
       % tic
        [f4,cn,E4,dE4,C4,H4] = ...
            invQnew(q0, B, Bu, Bv, res, stepsize, cutoff, itermax, cn);       
           
       % toc
        clear B;
        clear Bu;
        clear Bv;        
    end 
    
    %% I am going to visualize how the surface has eveolved over the optimization iterations
    f4_multires{i} = rotate3D(f4, O');
    cn_multires{i} = cn;
    
    E{i} = E4;
    if nargout == 4,
        dE{i} = dE4;
    end
end

f4_multires0 = cell(nres, 1);
for i=1:nres,
    f4_multires0{i} = permute(reshape( f4_multires{i}, [d, params.RES(i, 1), params.RES(i, 2)]), [2, 3, 1]);
end
f4_multires = f4_multires0;


%% visualization