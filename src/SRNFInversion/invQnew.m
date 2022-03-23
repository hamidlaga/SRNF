function [f,c,E,dE,C,H] = invQnew(q,B,Bu,Bv,res,stepsize,cutoff,itermax,c_init)

[Theta,Phi] = genGridSphr(res);
mysmall = 0;
sinPhi  = genSinPhi(Theta,mysmall);
dtheta  = Theta(1, 2) - Theta(1, 1);
dphi    = Phi(2, 1)   - Phi(1, 1);

dA = dtheta*dphi;

d  = 3; 
N  = prod(res);
nB = size(B,3);

%initialize with a sphere
C = []; %zeros(nB,itermax);
if (nargin == 8) || (numel(c_init)== 0)
    c = zeros(nB,1);
    c(1) = 1;
    c(5) = 1;
    c(9) = 1; 
elseif (nargin == 9)
    c = c_init;
end
% C(:,1) = c;
f = genSurfR3(c, B);
E    = zeros(1,itermax);
iter = 1;

[fu,fv]      = partialF(f,sinPhi,d,N,res);
[Qf,nf,~,r1] = map2qnew(fu,fv,N);    

%%
w        = Qf - q;
intgrand = innerS2(w,w,sinPhi,N); 
E(iter)  = intgrand*dA;
fprintf('Iteration %d with energy %f \n',iter, E(iter));

% h2 = figure(1200); clf;
% dispSurfR3(f,res,3);
% pause

dE = zeros(1,itermax);

for iter = 2:itermax
    
    gradC = gradInvQnew(Bu, Bv, fu, fv, nf,w,r1,sinPhi,dA, N, nB);   
    
    
    % update the surface
    c_n = c - stepsize.*gradC(:);
    f_n = genSurfR3(c_n,B);
    [fu_n,fv_n]      = partialF(f_n,sinPhi,d,N,res);
    [Qf,nf_n,~,r1_n] = map2qnew(fu_n,fv_n,N);
    
    w_n = Qf - q;   % the error ??
   
    
    
    
    intgrand = innerS2(w_n,w_n,sinPhi,N);
    E(iter)  = intgrand*dA;    
    dE(iter) = gradC' * gradC;

   
    if (E(iter) > E(iter-1)), %  && ~changed,
        stepsize = stepsize/2;         
        E(iter)  = E(iter-1);
        % disp('continuing');
        continue;
     else
        c  = c_n;
        f  = f_n;
        fu = fu_n;
        fv = fv_n;
        
        nf = nf_n;
        r1 = r1_n;
        w  = w_n;
    end
    
   %  C(:,iter) = c;    
    if (dE(iter) < cutoff), %  || (abs((E(iter) - E(iter-1)) / E(iter-1)) < cutoff)       
        break
    end 
    
    if abs(E(iter) - E(iter-1))/E(iter) < cutoff,
        break;
    end
    
    
end
E  = E(1:iter);
dE = dE(2:iter-1);
%C  = C(:,1:iter);

H  = [];

end
