function Mn = normalize_translation(M, Theta)
  
[n, m, d]   = size(M); 
if nargin < 2,
    [Theta, ~] = genGridSphr([n, m]);
end

[A1,multfact1,~] = area_surf_closed(M);
Mn = center_surface_closed(M, multfact1,A1,Theta);