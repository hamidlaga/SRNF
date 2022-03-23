%
% w : 3 x N
function [gradC] = gradInverseQnew_Hamid(Bu, Bv, ... 
                                   fu, fv, ...
                                   nf, w, r1, ...
                                   sinPhi, dA)
                               
% double *wdnf, *r5, nb[3], intgrand=0.0, wdnb=0.0, ndnb=0.0;
%     
%     wdnf = new double[N];
%     r5   = new double[N];
    
%disp('Here....');

MYLARGE = 10000.0;
nB      = size(Bu, 3);
N       = size(Bu, 2);

wdnf = sum(w .* nf, 1);     % <w, nf>
r5   = r1 .^ 5;
r5 (r5 > MYLARGE) = 0;

gradC = zeros(nB, 1);
for i=1:nB,

    intgrand = 0.0;
    
    bu = squeeze(Bu(:,:, i));
    bv = squeeze(Bv(:,:, i));
    
    for j=1:N,
       nb = crossDiff( fu(:, j), fv(:, j),  bu(:, j), bv(:, j));
       
       wdnb = dot(w(:, j), nb);
       ndnb = dot(nf(:, j), nb);
       
       intgrand = intgrand + (2*wdnb * r1(j) - ndnb * wdnf(j) * r5(j)) * sinPhi(j);
    end
    
    gradC(i) = gradC(i) + intgrand * dA;
    
end


%     for (int j=0; j<N; j++) {
%         wdnf[j] = dot3(w+j*3,nf+j*3);
% 
%         
%         r5[j] = pow(r1[j], 5);
%         if (r5[j] > MYLARGE)
%             r5[j] = 0.0;
%     }
%     
%     for (int i=0; i<nB; i++) {   
%         intgrand = 0.0;
%         for (int j=0; j<N; j++) {
% 			
%             crossDiff(nb,fu+j*3, fv+j*3, Bu+i*3*N+j*3, Bv+i*3*N+j*3);
%             
%             wdnb = dot3(w+j*3,nb);
%             ndnb = dot3(nf+j*3,nb);
% 			
%             intgrand += (2*wdnb*r1[j] - ndnb*wdnf[j]*r5[j])*sinPhi[j];
%         }
%         gradC[i] = intgrand*dA;
%     }
%     delete [] wdnf;
%     delete [] r5;
% }

function z = crossDiff(fu, fv, gu, gv) 
z(1) = fu(2)*gv(3)-fu(3)*gv(2) + gu(2)*fv(3)-gu(3)*fv(3);
z(2) = -fu(1)*gv(3)+fu(3)*gv(1) -gu(1)*fv(3)+gu(3)*fv(2);
z(3) = fu(1)*gv(2)-fu(2)*gv(1) + gu(1)*fv(2)-gu(2)*fv(1);


