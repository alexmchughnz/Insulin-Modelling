FolderName = "C:\Users\adm181\Google Drive\Work\PhD\Insulin Modelling\4. Detemir Resimulation\Plots";   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  set(0, 'CurrentFigure', FigHandle)
  set(gcf, 'Position', [0 0 800 600])
  
  AxisHandle = FigHandle.Children(end);
  FigName   = get(AxisHandle.Title, 'String');
  FigName = matlab.lang.makeValidName(FigName);
  savefig(FigHandle, fullfile(FolderName, 'fig', FigName + ".fig"));
  saveas(FigHandle, fullfile(FolderName, 'png', FigName + ".png"));
  
  pdfFile = fullfile(FolderName, 'pdf', FigName + ".pdf");
  saveas(FigHandle, pdfFile);
%   system(sprintf("pdfcrop %s %s", pdfFile, pdfFile));
end