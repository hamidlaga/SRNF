function multiRes_n = normalize_translation_multiresSurf(multiResM)

multiRes_n = multiResM;

for i=1:size(multiResM, 1),
    M = squeeze(multiResM(i, :,:,:));
    
    multiRes_n(i, :,:,:) = normalize_translation(M);
    
end


