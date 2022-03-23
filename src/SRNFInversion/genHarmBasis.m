function [B,Bu,Bv,Cij] = genHarmBasis(res,L,flag_Cij)

if nargin < 3
    flag_Cij = 0;
end

[Theta,Phi] = genGridSphr(res, 0.0);
dtheta = Theta(1,2)- Theta(1,1);
dphi   = Phi(2, 1) - Phi(1, 1);

% N = prod(res);
% sinTheta = genSinPhi(Theta, 0.00);


[~,~,~, Y] = SHarmBasis(res, L);
Y = Y(:,:, 2:end);


Y = GramSchmidt_Hamid(Y, res);
N = size(Y,3); 

B = SHarmBasisR3(reshape(Y, [prod(res), N]));

%% normalize the norm of the basis (unit surface area
% for k=1:size(B, 3),
% 
%     B1 = squeeze(B(:,: , k));
%     A2 = innerS2(B1, B1, sinTheta, N);% * dtheta * dphi;
%     if A2>0.00001,
%         B(:,: , k) = B(:,: , k) / sqrt(A2);
%     end
% end

if nargout > 1,
    [Bu,Bv] = BasisDiffR3(B,res);
end

if (flag_Cij)

    disp('(Basis) Compute cross products and parameters...')
    Cij = compBasisCross(Bu,Bv);

    clear Y b nY N Bu Bv
else
    Cij = [];
end
