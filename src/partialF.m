function [Fu,Fv] = partialF(f,sinTheta,numdim,npt,res)

% dtheta = 2*pi/(res(1)-1);
% dphi   =   pi/(res(2)-1);

dphi   = 2*pi/res(1); % /(res(1)-1);
dtheta =   pi/res(2); % (res(2)-1);


Fu = zeros([numdim,npt]);
Fv = zeros([numdim,npt]);
F  = reshape(permute(f,[2,1]),[res,numdim]);


for nd=1:numdim
    fu = DGradient(F(:,:,nd),dphi,1);  % DGradient(F(:,:,nd),dtheta,2);
    
    fv = DGradient(F(:,:,nd),dtheta,2);% DGradient(F(:,:,nd),dphi,1);

    
   %[fv, fu] = gradient(F(:,:,nd),dtheta,dphi); 
   Fu(nd,:) = fu(:);
   Fv(nd,:) = fv(:)./sinTheta(:);
end

