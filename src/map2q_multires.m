function Q_multires = map2q_multires(M_multires,  RES, mysmall )

d    = 3;                           % dimension of the data
nres = size(RES, 1);
Q_multires = cell(1, nres);

for i=1:nres,
    res = RES(i, 1:2);    
    [Theta, ~] = genGridSphr(res, mysmall);    
    sinTheta   = genSinPhi(Theta, mysmall);
    N       = prod(res);    
    
    [Mu, Mv] = partialF(M_multires{i}, sinTheta, d, N, res);
    [Q_multires{i}, ~, ~, ~] = map2qnew(Mu, Mv, N);    
end