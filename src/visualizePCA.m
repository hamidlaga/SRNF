function visualizePCA(Mu, B, res)
%
% Visualizes the mean and the modes of variation
%

%% The mean
M = reshape(Mu, [3, res]);
dispSurfR3(mySurf(M), res, 3);

n = size(B, 3);

eigenVals = ones(1, n);

for modeId = 1:n,  % Let's visualize 10 modes

    % close all;
    eVal  = sqrt(eigenVals(modeId));        % standard deviation

    f = -2:0.5:2; % [-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3];
    Variations =  f* eVal;

    h = figure(modeId * 1000);
    clf;
    hold on;
    colors = {'.r' ; '.g' ; '.b'; 'xg'; 'xr'};
    
    shift_x = 1;
    for i=1:length(Variations)
        sc = Variations(i);
        
        S  = reshape(Mu + sc * squeeze(B(:,:,  modeId)), [3, res] ); 
        
        f0 = reshape(S, [3,res]);
        S  = permute(f0,[2,3,1]);

        surface(squeeze(S(:, :, 1)) + shift_x * i, squeeze(S(:, :, 2)), squeeze(S(:, :, 3)));
    end   
  
    axis equal;
    cameramenu;
    axis off;
    
    lighting phong;
    light('Position',[0 0 1],'Style','local')
%    camlight('headlight');
    shading 'flat';
    set(gcf, 'Renderer', 'zbuffer');
    
	% savefig(h, [outDir 'Hamate_mode_' num2str(modeId) '.fig']);
    disp('Press any key to visualize the next mode of variation ...');
    pause; % (.1)
    
end