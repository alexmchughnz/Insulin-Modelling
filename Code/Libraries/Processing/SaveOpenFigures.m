function SaveOpenFigures(Trial, trialPath)


FolderName = Trial.Config.PLOTPATH;   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    % Find
    FigHandle = FigList(iFig);
    set(0, 'CurrentFigure', FigHandle)
    FigHandle.Units = 'pixels';
    FigHandle.Position = [0 0 600 400];
    
    AxisHandle = FigHandle.Children(end);    
    
    % Edit    
    FigName = MakeValidName(FigHandle.Name);
    
    if isfield(AxisHandle, 'Legend') && ~isempty(AxisHandle.Legend)
        legend('Location','southoutside')
    end
    
    % Save    
    figDir = fullfile(FolderName, trialPath, "fig");
    if ~isfolder(figDir)
        mkdir(figDir);
    end
    savefig(FigHandle, fullfile(figDir, FigName + ".fig"));
    
    pngDir = fullfile(FolderName, trialPath, "png");
    if ~isfolder(pngDir)
        mkdir(pngDir);
    end
    saveas(FigHandle, fullfile(pngDir, FigName + ".png"));
end
end