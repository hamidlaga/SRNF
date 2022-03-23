function [theta,phi,r]=cartesian_to_sph_m(x,y,z)

r    =sqrt(x.^2+y.^2+z.^2);
theta=acos(z./sqrt(x.^2+y.^2+z.^2));
phi  =atan2(x,y);

[I1,I2]=find(phi<0);

for j=1:length(I1)
    phi(I1(j),I2(j))=phi(I1(j),I2(j))+2*pi;
end
