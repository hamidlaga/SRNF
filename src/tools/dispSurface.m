function dispSurface(F, mycam)

if nargin <2
    mycam = [0, 0, 30];
end

[h, w, ~] = size(F);
cmap = hsv;
cmap = [cmap; 1 0.95 0.9; 0.7 0.7 0.7];
ncolors = size(cmap, 1);

surf(squeeze(F(:,:,1)), squeeze(F(:,:,2)), squeeze(F(:,:,3)),  ...
     ncolors * ones( h, w), 'SpecularExponent',100, ...
    'FaceAlpha', 1); % .5);
cameramenu;

shading flat
colormap(cmap);

% lighting flat 
axis equal off;
light; lighting phong; 

%view(0,90);
% campos(mycam);
% camup([0     0    -1]);