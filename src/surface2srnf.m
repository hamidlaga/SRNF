function [Q, Qu, Qv] = surface2srnf(S, Theta)
%
% It returns Q, the SRNF
% If nargout > 1, it also returns Qu, the SRVF along the u, and Qv, the SRVF along the v directions
%
% u: 0 to 2PI
% v: 0 to PI
%%
    
[resolution(1), resolution(2), ~] = size(S);

myS     = mySurf(S, 0);
mysmall = 0.01;
if nargin < 2
    [Theta, ~] = genGridSphr(resolution);
end

sinPhi  = genSinPhi(Theta, mysmall);
N       = prod(resolution);
d       = 3;

[Mu, Mv] = partialF(myS, sinPhi, d, N, resolution);
[myQ,  ~, ~, ~] = map2qnew(Mu, Mv, N);

Q  = permute(reshape(myQ, [3, resolution]), [2,3,1]);
Qu = zeros(size(Q));
Qv = zeros(size(Q));

if nargout > 1   
    
    %% Closed curve
    for i=1:size(Qu, 2)
        f = squeeze(S(:, i, :)); 
        % f = f';
        % f = f - repmat(f(:, 1), 1, size(f, 2)); % not needed
        q = curve_to_srvf_closed(f', 0); %, 2*pi); % q = curve_to_q(f', 0); %
        Qu(:, i, :) = q(1:3, :)';
                
%          curve = srvf_to_curve_closed(q);
%          
%          figure(101); clf;
%          plot3(f(:, 1), f(:, 2), f(:, 3), '-r'); hold on;
%           plot3(curve(1, :), curve(2, :), curve(3, :), '-b'); hold on;
%           pause
         
            
    end

    %% Open curve
    for j=1:size(Qv, 1)
        f = squeeze(S(j, :, :)); 
        %f = f';
        % f = f - repmat(f(:, 1), 1, size(f, 2));  % not needed
        q = curve_to_srvf(f', 0); %, pi); % q = curve_to_q(f', 0); 
        Qv(j, :, :) = q(1:3, :)';               
        
%          curve = srvf_to_curve(q);
%          
%          figure(101); clf;
%          plot3(f(:, 1), f(:, 2), f(:, 3), '-xr'); hold on;
%           plot3(curve(1, :), curve(2, :), curve(3, :), '-ob'); hold on;
%           pause


    end
end