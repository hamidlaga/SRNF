function Sn = resampleSphericalGrid(S0, N, toflip)
%
% samplingMode: 1 --> interpolation, 2 --> each two points pick 1
% 
% This function also replaces the NaN values (if any) by the average of the
% neighoring values

if nargin < 3,
    toflip = 0;
end

%% Clearning the NaNs
S0 = cleanNaNSurface(S0);

% make the data in same convention as Sebastian
if toflip
    Stmp = permute(S0,[2,1,3]);
else
    Stmp = S0;
end

fn = MakeClosedGrid(Stmp, N);       % 

% back to my convention
Sn = [];
if toflip   
    Sn = permute(fn, [3 1 2]);
else
   Sn = fn; 
end