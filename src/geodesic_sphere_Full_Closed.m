% This function calculates a geodesic on the sphere between q1 and q2
%
% range is the sampling of the cuve, which can be 1 (if t\in[0, 1]), 2pi if
% t\in[0, 2\pi],, etc
%
% set isClosed to 1 for closed curves
%
function X = geodesic_sphere_Full_Closed(q1,q2,stp,isClosed)

if nargin<4
    isClosed = 0;
end

if nargin < 3
    stp = linspace(0, 1, 7);
end

n = size(q1, 2);
X = zeros(3, n, numel(stp));

%% Global
norm1 = sqrt(InnerProd_Q_closed(q1, q1 )); %,range));
norm2 = sqrt(InnerProd_Q_closed(q2, q2 )); %,range));

q1n   = q1 / norm1;
q2n   = q2 / norm2;
theta = acos(InnerProd_Q_closed(q1n,  q2n ));%,  range));

if theta < 0.0001
    for t=1:stp+1
       X(:,:, t) = q1;
    end    
    return;
end

%% Individual vectors
% Norm1 = sqrt(sum(q1.^2, 1 ));       % these are the norms of individual vectors
% Norm2 = sqrt(sum(q2.^2, 1 ));
% 
% q1nn  = q1 ./ repmat(Norm1, 3, 1);
% q2nn  = q2 ./ repmat(Norm2, 3, 1);
% 
% Theta = acos(sum(q1nn .* q2nn, 1));
% Theta = repmat(Theta, 3, 1);
% sinTheta = sin(Theta);

X(:,:, 1) = q1;
X(:,:, numel(stp)) = q2;

for t=2:numel(stp)-1
    tau = stp(t); % (t-1)/stp; 
%    X(:,:, t) = (1-tau).*q1 + tau*q2;      % Linear

    %% use either this one - this one is better but still has some artifacts
%     X(:,:, t) =  (sin((1-tau)*Theta).*q1nn + sin(tau*Theta).*q2nn) ./ sinTheta; % Recover the orientation
%     X(:,:, t) = ((1-tau)*Norm1 + tau*Norm2) .* X(:,:, t);
     
    %% Or this one - check the difference - this one deforms the entire surface as one vector
    %  Thus theoretically it should not be affected too much if two specific normals
    %  have an angle larger than 180 degrees, since it will find the optimal overall direction.
    X(:,:, t) = (sin((1-tau)*theta)*q1n + sin(tau*theta)*q2n)/sin(theta); 
    
    if isClosed==1  
        %% this makes sure that the curve is closed
        X(:,:,t) = myProjectC(X(:,:,t));
    end  
    X(:, :, t) = ((1-tau)*norm1 + tau*norm2) * X(:, :, t);
    % (sin((1-tau)*theta)*norm1 + sin(tau*theta)*norm2)/sin(theta)  * X(:,:, t);  ;  %
end

return;









