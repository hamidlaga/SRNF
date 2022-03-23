function M_multires = prepareMultiResSurface(multiResM, params)

% Align the surface to the sphere
multiResM = normalize_translation_multiresSurf(multiResM);
multiResM = normalize_scale_multiresSurf(multiResM);

% Register it to a spehre
[~, O] = alignSurfaceGridWithSphere(squeeze(multiResM(end, :,:,:)), params.surfaceCodeDir, params.currDir); 
for i=1:size(multiResM, 1),
    multiResM(i, :,:,:) = rotateGrid(squeeze(multiResM(i, :,:,:)), O);
end

% Resample the multi-resolution representation to the desired resolutions
M_multires = cell(1, params.nres);
for i=1:params.nres,
    res = params.RES(i, 1:2);    
    M_multires{i} = resampleSphericalGrid(squeeze(multiResM(i, :,:,:)), res(1), 0);       
end