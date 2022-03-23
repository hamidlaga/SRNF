function p = srvf_to_curve(q)
%
qn    = q(:, 1:end);
[n,T] = size(qn);

qnorm =  sqrt(sum(qn.^2, 1));

for i = 1:n
    p(i,:) = [cumtrapz( qn(i,:).*qnorm )/(T) ] ;
end

