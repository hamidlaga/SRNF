%
% q2 - q1 is teh deformation to transport to q3
% qn = q3 + s* (q2 - q1) but I need to take into account the fact that the
% space of SRNFs is not linear
%
function qn = transportSRNFdeformation(q1, q2, s, q3)

q1 = mySurf(q1, 0);
q2 = mySurf(q2, 0);
q3 = mySurf(q3, 0);

qn = q1*(1-s) + q2*s;

[d, n] = size(q1);
vold   = [];

m = sqrt(n);

for i=1:n,
   v1 = q1(:, i);
   v2 = q2(:, i);
   v3 = q3(:, i);
%    
% disp(v1)
% pause
% disp(v2)
% pause
   
   a1 = dot(v1, v1); % local area
   a2 = dot(v2, v2);
   a3 = dot(v3, v3);
   
   an = a3 + s * (a2 - a1); % a3 * s * (a2 / a1);
   
   v1 = v1 / sqrt(a1);
   v2 = v2 / sqrt(a2);
   v3 = v3 / sqrt(a3);
    
   % theta = acos(dot(v1, v2));   % angle between v1 and v2
   
   
   r = vrrotvec(v1,v2);
   r(4) = s * r(4);
   O = vrrotvec2mat(r);
   
%    cth = dot(v1, v2);       % cos
%    sth = norm(cross(v1, v2));     % sin
%    theta = atan2(sth, cth);
%    thetan = s*theta;
%    cth = cos(thetan);
%    sth = sin(thetan);
%    
%    O   = [cth, -sth; ...
%           sth, cth];

   vn = O * v3(:);  
   
   vn = v3 + s * (v2 - v1);
   vn = vn / norm(vn);
   
%    theta_n = s * theta;         % the rotation to apply to v3
%    
%    vn = (1 - s) * v1 / norm(v1) + s * v2/norm(v2);
%    
%    vn = vn/norm(vn);
   

   qn(:, i) = sqrt(an) * vn;

  
end

 qn = permute(reshape( qn, [d, m, m]), [2, 3, 1]);
