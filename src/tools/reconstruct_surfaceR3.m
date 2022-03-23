function [c_new,f_new] = reconstruct_surfaceR3(f_orig,myInner,B,res, Mu)

if nargin < 5,
    Mu = zeros(3, prod(res));
end

% disp('(Surface) Reconstructing surface...')
myView = [53,18];

[d,nPt,nB]  = size(B);
[Theta,Phi] = genGridSphr(res);

f_new = zeros(size(f_orig));
c_new = zeros(nB, 1);

d = 3;
dtheta = Theta(1, 2) - Theta(1, 1); % 2*pi/(res(1)-1);
dphi   = Phi(2, 1)   - Phi(1, 1);   %   pi/(res(2)-1);

for j=1:nB
    tmp      = squeeze(sum( (f_orig - Mu) .* squeeze(B(:,:,j)),1));
    intgrand = reshape(tmp(:),res).* sin(Theta);
    c_new(j) = sum(intgrand(:)) * dtheta*dphi;
end
% disp(c_new)
% pause
f_new = genSurfR3(c_new, B) + Mu;

% show_fig = 1;
% if (show_fig)
% 	figure;    
%     %for i=1:d
%         subplot(2, 1, 1); dispSurfR3(f_orig, res, d); %dispSurf(Theta,Phi,reshape(f_orig(i,:)',res),1);
%         subplot(2, 1, 2); dispSurfR3(f_new, res, d); %dispSurf(Theta,Phi,reshape(f_new(i,:)',res),1);
%     %end
%     pause
% end

show_fig = 0;
if (show_fig)
    figure; clf
    subplot(1,2,1);
    dispSurfR3(f_orig,res,d,myView);
    subplot(1,2,2);
    dispSurfR3(f_new,res,d,myView);
    pause
end


end