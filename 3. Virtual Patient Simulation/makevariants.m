SIValues = [2.5e-4, 5e-4, 10.8e-4, 18e-4];

variants = cell(1,length(SIValues));
for ii = 1 : length(variants)
    V.SI = SIValues(ii);
    V.results = struct();
    
    variants{ii} = V;
end

save('variants.mat', 'variants')
disp('Variants created.')
clear