function FigHandle = PanelFigures(numFigures, monitor)

persistent n;
if (isempty(n))
    n = 1;
end

G = num2cell(get(groot, 'MonitorPositions'));

if (isempty(monitor) || (size(G,1) < monitor))
    monitor = 1;
end

screensize = num2cell(get(groot, 'MonitorPositions'));
[x0, y0, w, h] = screensize{monitor, :};
FigHandle = figure(n);
FigHandle.Position = [x0 + (n-1)/numFigures*w, y0, w/numFigures, 0.9*h];
n = n + 1;

end

