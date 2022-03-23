function p = srvf_to_curve_general(q, X)
%
% X is the sampling
%
qn    = q(:, 1:end);
[n,T] = size(qn);
qnorm = sqrt(sum(qn.^2, 1));

for i = 1:n
    p(i,:) = cumtrapz(X,  qn(i,:).*qnorm )/(X(end) - X(1)); %(T) ] ;
end