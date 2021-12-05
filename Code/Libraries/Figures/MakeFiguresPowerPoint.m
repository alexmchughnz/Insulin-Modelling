function MakeFiguresPowerPoint(Trial)
% Save all patient plots in a powerpoint.


if ~ispc
    disp("Can only make PowerPoint on PC.")
    return
end


numPatients = numel(Trial.patientSet);

figDir = fullfile(Trial.Config.PLOTPATH, Trial.outputPath);
if ~isfolder(figDir)
    mkdir(figDir);
end

recipeName = split(Trial.outputPath, "\");
recipeName = recipeName{end};

pptDir = fullfile(figDir, Trial.source+"_"+recipeName);


for pp = 1:numPatients
    figsPerSlide = 4;
    
    P = Trial.patientSet{pp};
    allFigs = P.figures;
    
    numSlides = ceil(numel(allFigs)/figsPerSlide);
    
%     saveppt2(pptDir, 'figure', 0, )
    for nn = 1:numSlides
        range = [1:figsPerSlide] + (nn-1)*figsPerSlide;
        range = range(range<=numel(allFigs));
    
        saveppt2(pptDir, 'figure', allFigs(range), 't', P.patientCode, 'scale', true, 'stretch', false)
    end
    
end
end

