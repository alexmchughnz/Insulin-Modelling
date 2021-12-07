function SaveFigures(Trial, P)
% Saves all figures belonging to a patient.


for ii = 1:length(P.figures)
    % Find
    F = P.figures(ii);
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
    figDir = fullfile(Trial.Config.PLOTPATH, Trial.outputPath, "fig");
    if ~isfolder(figDir)
        mkdir(figDir);
    end
    savefig(F, fullfile(figDir, figName + ".fig"));
    
    pngDir = fullfile(Trial.Config.PLOTPATH, Trial.outputPath, "png");
    if ~isfolder(pngDir)
        mkdir(pngDir);
    end
    saveas(F, fullfile(pngDir, figName + ".png"));
    
    % Close
    if Trial.Config.CLOSEALLFIGURES
        close(F)
    end
end
end