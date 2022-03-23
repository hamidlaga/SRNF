function M = invertSRNF_Analytical(Q)
%
% This function is only suitable for star-shaped surfaces
% (e.g., the surfaces in ./sample_surfaces/Lthalamus.mat)
%

%% 1. Compute the E
[res(1), res(2), ~] = size(Q);

% 1.1. Get a sphericla grod
[Theta, Phi] = getSphericalGrid(res(1), 0.01);
[X, Y, Z]    = spherical_to_cart_m(Theta, Phi);
E(:,:, 1) = X;
E(:,:, 2) = Y;
E(:,:, 3) = Z;

%% 2. Compute Q^r
Qr = sum(Q .* E, 3);

%% 3. Reconstruct the surface
normQ =  sqrt(sum(Q.*Q, 3));
M     = repmat(sqrt(normQ .* abs(Qr)), [1, 1, 3]) .* E;

for i=1:size(M, 1)
    for j=1:size(M, 2)
        p = squeeze(M(i, j, :));
%         if norm(p) > 2,
%             M(i, j, :)=[0,0,0];
%         end
        
        if Qr(i, j)<0,
            M(i, j, :)=[0,0,0];
        end
    end
end

