% This function calculates a geodesic on the sphere between q1 and q2
%
% range is the sampling of the cuve, which can be 1 (if t\in[0, 1]), 2pi if
% t\in[0, 2\pi],, etc
% 
% Computes geodesic between (q1, n1) and (q2, n2)
% q is the SRVF
% n is the SRNF thus they are orthogonal
%
% Takes cross q and n which results in another orthogonal vector p. Once
% normalised (q, n, p) is a rotation matrix - thus the geodesic becomes a
% geodesci between rotation matrices
%
%
function X = geodesic_rots_Closed(q1,n1, q2, n2, Tau, toProject)

if nargin<6
    toProject = 0;
end

if nargin < 5
    Tau = linspace(0, 1, 7);
end


[~,n] = size(q1);

p1 = zeros(3, n);
p2 = zeros(3, n);

for i=1:n       
    p1(:, i) = cross(q1(:, i), n1(:, i));        
    p2(:, i) = cross(q2(:, i), n2(:, i));
end

norm1 = sqrt(InnerProd_Q_closed(q1, q1)); %,range));
norm2 = sqrt(InnerProd_Q_closed(q2, q2)); %,range));

q1n = q1 / norm1;
q2n = q2 / norm2;
theta = acos(InnerProd_Q_closed(q1n,  q2n)); %,  range));


Norm1 = sqrt(sum(q1.^2, 1 ));
Norm2 = sqrt(sum(q2.^2, 1 ));

% Theta = acos(sum((q1 ./ repmat(Norm1, 3, 1)) .*  (q2 ./ repmat(Norm2, 3, 1) ), 1) );
% sinTheta = sin(Theta);

% Norm_n1 = sqrt(sum(n1.^2, 1 ));
% Norm_n2 = sqrt(sum(n2.^2, 1 ));

if theta > 0.0001
    
    for i=1:n            
        v1 = squeeze(q1(:, i)); v1 = v1 / norm(v1);
        v2 = squeeze(n1(:, i)); v2 = v2 / norm(v2);
        R1 = [v1, v2, cross(v1, v2)];
        
        v1 = squeeze(q2(:, i)); v1 = v1 / norm(v1);
        v2 = squeeze(n2(:, i)); v2 = v2 / norm(v2);
        R2 = [v1, v2, cross(v1, v2)];
        
        Os = interpolateRots(R1, R2, Tau);       % interpolation

        for t=1:numel(Tau)
            X(:, i, t) = ((1 - Tau(t))*Norm1(i) + Tau(t)*Norm2(i))*Os{t}(:, 1) ; 
        end
    end
%     
    if toProject==1
        for t=1:numel(Tau)
        
            %% this makes sure that the curve is closed
            X(:,:,t) = myProjectC(X(:,:,t));
        end
    end
else
    for t=1:stp+1
        X(:,:, t) = q1;
    end
end

return;
