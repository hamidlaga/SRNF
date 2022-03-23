function [x,y,z] = dispSurf(f,Theta,Phi,is_sphr)

if nargin > 3 && is_sphr
    [x,y,z]=sphr2cart(Theta,Phi,squeeze(f));
else
    x = squeeze(f);
    y = squeeze(Theta);
    z = squeeze(Phi);
end

axis equal off; hold on;
surf(x,y,z); 
light; lighting phong; 
% camzoom(1.3);
view([30 25]);

% cmap = hsv( prod(size(x)));
% colormap([cmap; 1 0.95 0.9]);
cameramenu

end