function Os = interpolateRots(R1, R2, st)


Q1 = qGetQ(R1);  % first quaternion
Q2 = qGetQ(R2);  % second quat

Os = cell(1, numel(st));
for i=1:numel(st)
   q =  slerp(Q1, Q2, st(i), .01);
   Os{i} = qGetR(q);
end

% Os{numel(st)} = R2;
%  id_Q = qGetQ(eye(3, 3));  % identity
%     Q = qGetQ(inv(R));        % quaternion from q1 to q2
%     Q_half = slerp(id_Q, Q, 0.5, .01);
%     R_half = qGetR(Q_half);