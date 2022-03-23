function [QQ, QQu, QQv] = computeGeodesicSRNF(Q1, Qu1, Qv1, Q2, Qu2, Qv2, st)
%
% It assumes that the surfqaces have been optimally reparameterized and registered. As such,
% this function does not find the optimal rotation / reparameterization
%
if nargin < 7
    st = linspace(0, 1, 3);          % the number of samples along the geodesic
end
n = numel(st);

res = size(Q1);

QQ  = zeros(n, res(1), res(2), res(3));
QQu = zeros(n, res(1), res(2), res(3));
QQv = zeros(n, res(1), res(2), res(3));

%% Geodesic between the Qs
for i=1:n,
    alpha = st(i);
    QQ(i, :, :, :) = (1- alpha) * Q1 + alpha * Q2;
end

%% geodesic between the Qus
for i=1:size(Qu1, 2)
   
    q1 = squeeze(Qu1(:, i, 1:3))';
    q2 = squeeze(Qu2(:, i, 1:3))';
    
    norm1 = sqrt(InnerProd_Q_closed(q1, q1)); %,range));
    norm2 = sqrt(InnerProd_Q_closed(q2, q2)); %,range));

    if (norm1 == 0) || (norm2 == 0)
        R = eye(3, 3);
    else        
        % R aligns q2n onto q1n
        [~,R] = Find_Best_Rotation(q1/norm1, q2/norm2);      
    end
        
    %R = eye(3, 3);
    Os = interpolateRots(eye(3, 3), inv(R), st);
    
    %% Use this
    % g = geodesic_sphere_Full_Closed(q1, R*q2, st,  1);    
    g = geodesic_sphere_Full_Closed_linear(q1, R*q2, st,  1);
    
    for j=2:n
        g(:,:, j) = Os{j}*squeeze(g(:,:, j));
    end
    QQu(:, :, i, :) = permute(g, [3 2 1]) ;
end

%% Geodesic between the Qvs
for i=1:size(Qv1, 1)
   
    q1 = squeeze(Qv1(i, :, 1:3))';
    q2 = squeeze(Qv2(i, :, 1:3))';   
    
    norm1 = sqrt(InnerProd_Q_closed(q1, q1)); %,range));
    norm2 = sqrt(InnerProd_Q_closed(q2, q2)); %,range));
    
    if (norm1 == 0) || (norm2 == 0)
        R = eye(3, 3);
    else 
        [~, R] = Find_Best_Rotation(q1/norm1, q2/norm2);    
    end
    %R = eye(3, 3);
    
    Os = interpolateRots(eye(3, 3), inv(R), st);
    
    %% Use this
    % g = geodesic_sphere_Full_Closed(q1, R*q2, st);
    g = geodesic_sphere_Full_Closed_linear(q1, R*q2, st);
    
    for j=2:n
        g(:,:, j) = Os{j}*squeeze( g(:,:, j));
    end
    QQv(:, i, :, :) = permute(g, [3 2 1]) ;
    
end

end