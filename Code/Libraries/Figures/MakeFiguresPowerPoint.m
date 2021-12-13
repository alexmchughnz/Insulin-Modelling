function MakeFiguresPowerPoint(Trial, patientSet)
% Save all patient plots in a powerpoint.


%% Setup
figsPerSlide = 4;

if ~ispc
    disp("Can only make PowerPoint on PC.")
    return
end

figDir = fullfile(Trial.Config.PLOTPATH, Trial.outputPath);
if ~isfolder(figDir)
    mkdir(figDir);
end

recipeName = split(Trial.outputPath, "\");
recipeName = recipeName{end};

pptDir = fullfile(figDir, Trial.source+"_"+recipeName+".ppt");

%% Save
for pp = 1:numel(patientSet)
    P = patientSet{pp};

    allFigs = P.figures;
    numSlides = ceil(numel(allFigs)/figsPerSlide);

    for nn = 1:numSlides
        range = [1:figsPerSlide] + (nn-1)*figsPerSlide;
        range = range(range<=numel(allFigs));

        title = sprintf("%s (%d/%d)", patientSet.patientCode, nn, numSlides);

        saveppt2(pptDir, 'figure', allFigs(range), 't', title, 'scale', true, 'stretch', false)
    end
end
end

