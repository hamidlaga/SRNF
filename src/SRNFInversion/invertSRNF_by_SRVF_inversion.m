function S = invertSRNF_by_SRVF_inversion(Q, Qu, Qv)

res = size(Q);

[Theta, Phi] = genGridSphr(res);

%% Let's invert the first vertical curve, i.e., Q(1, :. :), which is an open curve
q = squeeze(Qv(1, :, :) )';
refCurve = srvf_to_curve_general(q, Theta(1, :)); 

%% Now let's invert the azimuth (lattitude) curves
S1 = zeros(res(1), res(2), 3);
for i = 1:size(Q, 2)   
    q  = squeeze(Qu(:, i, :));
    q  = q';
    curve = srvf_to_curve_closed(q); % \general(beta', Phi(:, 1));    
    curve = curve + repmat(refCurve(:, i) + curve(:, 1), 1, size(curve, 2));

    S1(:, i, :) = reshape(curve', [1, res(1), 3]);  
end
S = S1;

% refCurve = srvf_to_curve(squeeze(Qu(:, 1, :) )'); 
% 
% %% Now let's invert the longitudinal curves
% S2 = zeros(res(1), res(2), 3);
% for i=1: size(Q, 1),   
%     beta = squeeze(Qv(i, 1:end, :));
%     curve = srvf_to_curve(beta'); %, refCurve(:, i)); 
%     curve = curve + repmat(refCurve(:, i)- curve(:, 1) , 1, size(curve, 2));
% 
%     S2(i, : ,:) = reshape(curve', [1 res(1), 3]); 
% end
% 
% S = (S1 + S2)/2;