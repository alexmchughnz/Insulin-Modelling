function FigHandle = PanelFigures(numFigures, monitor)

persistent n;
if (isempty(n))
    n = 1;
end

G = num2cell(get(groot, 'MonitorPositions'));

if (isempty(monitor) || (size(G,1) < monitor))
    n = 1;
end

screensize = num2cell(get(groot, 'MonitorPositions'));
[x0, y0, w, h] = screensize{3, :};
FigHandle = figure(gcf().Number + 1);
FigHandle.Position = [x0 + (n-1)/numFigures*w, y0, w/numFigures, 0.9*h];
n = n + 1;

end

