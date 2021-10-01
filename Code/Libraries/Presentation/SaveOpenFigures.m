function SaveOpenFigures(T, trialPath)

global CONFIG

FolderName = CONFIG.PLOTPATH;   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    % Find
    FigHandle = FigList(iFig);
    set(0, 'CurrentFigure', FigHandle)
    FigHandle.Units = 'pixels';
    FigHandle.Position = [0 0 600 400];
    
    AxisHandle = FigHandle.Children(end);    
    
    % Edit    
    FigName = MakeValidName(T.source + FigHandle.Name);
    
    if isfield(AxisHandle, 'Legend') && ~isempty(AxisHandle.Legend)
        legend('Location','southoutside')
    end
    
    % Save    
    figDir = fullfile(FolderName, 'fig', trialPath);
    if ~isfolder(figDir)
        mkdir(figDir);
    end
    savefig(FigHandle, fullfile(figDir, FigName + ".fig"));
    
    pngDir = fullfile(FolderName, 'png', trialPath);
    if ~isfolder(pngDir)
        mkdir(pngDir);
    end
    saveas(FigHandle, fullfile(pngDir, FigName + ".png"));
end
end