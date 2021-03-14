function [] = SaveOpenFigures(subfolder)

global CONFIG

if isempty(subfolder)
    subfolder = "";
end

FolderName = CONFIG.PLOTPATH;   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    % Find
    FigHandle = FigList(iFig);
    set(0, 'CurrentFigure', FigHandle)
    FigHandle.Units = 'pixels';
    FigHandle.Position = [0 0 400 300];
    
    AxisHandle = FigHandle.Children(end);
    
    
    % Edit    
    FigName   = get(AxisHandle.Title, 'String');
    FigName = matlab.lang.makeValidName(FigName);        
    
    if (CONFIG.PUBLISHPLOTS)
        title('')
    end
    
    if ~isempty(AxisHandle.Legend)
        legend('Location','southoutside')
    end
    
    % Save
    
    figDir = fullfile(FolderName, 'fig', subfolder);
    if ~isfolder(figDir)
        mkdir(figDir);
    end
    savefig(FigHandle, fullfile(figDir, FigName + ".fig"));
    
    pngDir = fullfile(FolderName, 'png', subfolder);
    if ~isfolder(pngDir)
        mkdir(pngDir);
    end
    saveas(FigHandle, fullfile(pngDir, FigName + ".png"));
    
    pdfDir = fullfile(FolderName, 'pdf', subfolder);
    if ~isfolder(pdfDir)
        mkdir(pdfDir);
    end
    saveas(FigHandle, fullfile(pdfDir, FigName + ".pdf"));
    %   system(sprintf("pdfcrop %s %s", pdfFile, pdfFile));
end
end