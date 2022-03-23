function dispSurfR3(f,res,d,myView)
%
% f is of size 3x npoints
%

% dim(f) = [d,n1*n2]
f0 = reshape(f, [d, res]);
F  = permute(f0,[2,3,1]);
[h, w, ~] = size(F);
cmap = hsv;
cmap = [cmap; 1 0.95 0.9; 0.7 0.7 0.7];
ncolors = size(cmap, 1);

surf(squeeze(F(:,:,1)), squeeze(F(:,:,2)), squeeze(F(:,:,3)),  ...
     ncolors * ones( h, w), 'SpecularExponent',100, ...
    'FaceAlpha', 1); % .5);
cameramenu;

% dispSurf(F(:,:,1),F(:,:,2),F(:,:,3));
% dispSurf([F(:,:,1),F(:,1,1)],[F(:,:,2),F(:,1,2)],[F(:,:,3),F(:,1,3)]); %

if (nargin < 4)
    % view([-6,24])
else
    view(myView)
end
shading flat
% colormap(cmap);

% lighting flat 
axis equal off;
light; lighting phong; 
