function sinPhi = genSinPhi(Phi,mysmall)

sinPhi = sin(Phi);
% mysmall = eps; %eps;
[n1,n2] = size(Phi);
for i1=1:n1
    for i2=1:n2
        if (sinPhi(i1,i2)>=0 && sinPhi(i1,i2) < mysmall)
            sinPhi(i1,i2) = mysmall;
        elseif (sinPhi(i1,i2)<0 && sinPhi(i1,i2) > -mysmall)
            sinPhi(i1,i2) = -mysmall;
        end
    end
end
sinPhi = sinPhi(:);
sinPhi = reshape(sinPhi,[1,n1*n2]);