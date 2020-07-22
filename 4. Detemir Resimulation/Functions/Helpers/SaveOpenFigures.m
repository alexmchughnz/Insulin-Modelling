FolderName = "C:\Users\adm181\Google Drive\Work\PhD\Insulin Modelling\4. Detemir Resimulation\Plots";   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  set(0, 'CurrentFigure', FigHandle)
  set(gcf, 'Position', [0 0 800 600])
  
  AxisHandle = get(FigHandle, 'CurrentAxes');
  FigName   = get(AxisHandle.Title, 'String');
  FigName = matlab.lang.makeValidName(FigName);
  savefig(FigHandle, fullfile(FolderName, FigName + ".fig"));
end