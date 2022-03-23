% This function calculates a geodesic on the sphere between q1 and q2
%
% range is the sampling of the cuve, which can be 1 (if t\in[0, 1]), 2pi if
% t\in[0, 2\pi],, etc
%
function X = geodesic_sphere_Full_Closed_linear(q1,q2, Tau, toProject)

if nargin<4
    toProject = 0;
end

norm1 = sqrt(InnerProd_Q_closed(q1, q1)); %,range));
norm2 = sqrt(InnerProd_Q_closed(q2, q2)); %,range));

q1n = q1 / norm1;
q2n = q2 / norm2;
theta = acos(InnerProd_Q_closed(q1n,  q2n)); %,  range));

% [~,n] = size(q1n);
% Norm1 = sqrt(sum(q1.^2, 1 ));
% Norm2 = sqrt(sum(q2.^2, 1 ));
% Theta = acos(sum((q1 ./ repmat(Norm1, 3, 1)) .*  (q2 ./ repmat(Norm2, 3, 1) ), 1) );
% sinTheta = sin(Theta);

if theta > 0.0001
    for t=1:numel(Tau)
        tau = Tau(t); 
        
        X(:,:, t) = (1 - tau).*q1n + tau.*q2n ; 
        if toProject==1
            %% this makes sure that the curve is closed
            X(:,:,t) = myProjectC(X(:,:,t));
        end
         X(:,:, t) = ((1-tau)*norm1 + tau*norm2) * X(:,:, t);
    end
else
    for t=1:stp+1
        X(:,:, t) = q1;
    end
end

return;
