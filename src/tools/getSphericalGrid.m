%
% Makes a spherical grid in theta and phi coordinates
% Theta: 0 to \pi
% Phi:   0 to 2\pi
%
function  [Theta, Phi] = getSphericalGrid(N, eps)

if nargin < 2,
    eps = 0;
end

theta = linspace(0+eps, pi-eps, N); % pi*[0:N]/N;
Theta = repmat(theta,N,1);

phi = linspace(0+eps, 2*pi-eps, N); % 2*pi*[0:N]/(N);
Phi = repmat(phi',1, N);
