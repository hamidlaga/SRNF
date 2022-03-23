function M = normalizeSurface(M, Theta)

[res(1), res(2), ~] = size(M);

if nargin < 2
    [Theta,~] = genGridSphr(res);
end
    
[A1,multfact1,~] = area_surf_closed(M);
M                = center_surface_closed(M, multfact1, A1, Theta);
[A1,~,~] = area_surf_closed(M);
M        = M /sqrt(A1);