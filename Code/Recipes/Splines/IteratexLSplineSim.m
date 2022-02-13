function PArray = IteratexLSplineSim(P)
% Recipe for fitting parameters to a model with a fixed xL and splines for
% nL. Values for xL are iterated.


%% Functions
xLArray = 0.5 : 0.05 : 0.9;

for ii = 1:numel(xLArray)
    xL = xLArray(ii);
    
    copyP = TagPatientCode(P, "xL = " +xL);
    
    copyP = SCLossSplineSim(copyP, {}, xL);
    
    PArray{ii} = copyP;
end

end

