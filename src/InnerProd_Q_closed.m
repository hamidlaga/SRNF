function val = InnerProd_Q_closed(q1, q2, range)

if nargin < 3,
    range = 1;
end

[n, T] = size(q1);
val    = trapz(linspace(0, range,T), sum(q1.*q2));

return;