clear
close all

load('../config', 'RESULTPATH');
resultsfile = @(filename) fullfile(RESULTPATH, filename);

files = {'grid nL[0 0.4]@0.02 xL[0.5 1]@0.02P1.mat', ...
    'grid nL[0 0.4]@0.02 xL[0.5 1]@0.02 bestP1.mat', ...
    'grid nL[0 0.4]@0.02 xL[0.5 1]@0.02 worstP1.mat'};

figure()
hold on
for ii = 1:length(files)
    F = files{ii};
load(resultsfile(F), ...
        'nLGrid', 'xLGrid', 'IResiduals');
    
delta = 1e-3;
    
xLRange = xLGrid(:, 1);
nLRange = nLGrid(1, :);

ppxL = griddedInterpolant(xLRange, IResiduals(:, 1));
ppnL = griddedInterpolant(nLRange, IResiduals(1, :));

xLTest = min(xLRange) : delta : max(xLRange);
nLTest = min(nLRange) : delta : max(nLRange);


[xLMinResidual(ii), xLIndex(ii)] = min(ppxL(xLTest));
[nLMinResidual(ii), nLIndex(ii)] = min(ppnL(nLTest));


plot([nLTest(1), nLTest(nLIndex(ii))], [xLTest(xLIndex(ii)), xLTest(1)])
xlabel("nL")
ylabel("xL")
end


[xi, yi] = polyxpoly();
