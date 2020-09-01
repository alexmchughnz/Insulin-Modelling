function [] = SaveOpenFigures(subfolder)
if isempty(subfolder)
    subfolder = "";
end

FolderName = "C:\Users\adm181\Google Drive\Work\PhD\Insulin Modelling\4. Detemir Resimulation\Plots";   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    set(0, 'CurrentFigure', FigHandle)
    %   set(gcf, 'Position', [0 0 800 600])
    
    AxisHandle = FigHandle.Children(end);
    FigName   = get(AxisHandle.Title, 'String');
    FigName = matlab.lang.makeValidName(FigName);
    
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