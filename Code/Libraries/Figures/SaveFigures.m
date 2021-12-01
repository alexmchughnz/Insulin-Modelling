function SaveFigures(Trial, P)
% Saves all figures belonging to a patient.


figureFolder = Trial.Config.PLOTPATH;

for ii = 1:length(P.figures)
    % Find
    F = P.figures{ii};
    figure(F)

    F.Units = 'pixels';
    F.Position = [0 0 600 400];
    
    ax = F.Children(end);    
    
    % Edit    
    figName = MakeValidName(F.Name);
    
    if isfield(ax, 'Legend') && ~isempty(ax.Legend)
        legend('Location','southoutside')
    end
    
    % Save    
    figDir = fullfile(figureFolder, Trial.outputPath, "fig");
    if ~isfolder(figDir)
        mkdir(figDir);
    end
    savefig(F, fullfile(figDir, figName + ".fig"));
    
    pngDir = fullfile(figureFolder, Trial.outputPath, "png");
    if ~isfolder(pngDir)
        mkdir(pngDir);
    end
    saveas(F, fullfile(pngDir, figName + ".png"));
end
end