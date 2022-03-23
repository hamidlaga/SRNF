function multiRes_n = normalize_scale_multiresSurf(multiResM)

multiRes_n = multiResM;


if iscell(multiResM)

    for i=1:size(multiResM, 1),

        [A1,~,~] = area_surf_closed(multiResM{i} );
        multiRes_n{i} = multiResM{i}  / sqrt(A1); 

    end
    
else    
    for i=1:size(multiResM, 1),

        [A1,~,~] = area_surf_closed(squeeze(multiResM(i, :,:,:)));
        multiRes_n(i, :,:,:) = multiResM(i, :,:,:) / sqrt(A1); 

    end

end