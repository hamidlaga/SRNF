function [q] = curve_to_srvf_closed(p, toNormalizeLength)

if nargin < 2
    toNormalizeLength = 1;
end

[n,N] = size(p);

for i = 1:n
    v(i,:) = gradient(p(i,:),1/N);
end

for i = 1:N
    L(i) = sqrt(norm(v(:,i),'fro'));
    if L(i) > 0.00001
        q(:,i) = v(:,i)/L(i);
    else
        q(:,i) = zeros(n,1);
    end
end

%q = q/sqrt(InnerProd_Q(q,q));

% %===
% [n,N] = size(p);
% for i = 1:n
%     v(i,:) = gradient(p(i,:),range/(N));
% end
% 
% % why dividing by N for the length????
% len = sum(sqrt(sum(v.*v)))/N;
% 
% %% Normalize for scale (if needed)
% if toNormalizeLength > 0
%     v = v/len;        % This makes the length of the curve equal 1
% end
% 
% %% Compute the SRVF
% q = zeros(size(p, 1), size(p, 2));
% 
% for i = 1:N
%     L = max(sqrt(norm(v(:,i))), 0.0001);
%     q(:,i) = v(:,i)/L;
% end