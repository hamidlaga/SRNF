function p = srvf_to_curve_closed(q)

[n,T] = size(q);
qnorm = zeros(1, T);
p     = zeros(n, T);

for i = 1:T
    qnorm(i) = norm(q(:,i),'fro');
end

X = linspace(0, 1, T);
for i = 1:n
    p(i, :) = cumtrapz(X,  q(i,:).*qnorm ); %/ T;
end


% p(:,end) = p(:,1);
return;