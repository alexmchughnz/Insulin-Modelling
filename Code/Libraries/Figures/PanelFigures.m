function PanelFigures(Trial, monitor)
% Iterate through all patients' plots and displays them in a grid.

if ~exist('monitor', 'var')
    monitor = 1;
end

CONST = LoadConstants();

numCols = numel(Trial.patientSet);
numRows = numel(Trial.patientSet{1}.figures);
for ii = 1:numCols
    P = Trial.patientSet{ii};
    handles = P.figures;

    for hh = 1:numRows
        F = figure(handles{hh});

%         % Grab patient info for figure.
%         ax = gca();
%         lgd = ax.Legend;
%         if ~isempty(lgd)
%         lgd.Location = 'best';
%         end

        % Get position indices.
        colIndex = ii;
        rowIndex = hh;

        % Position figure.
        F.Units = 'pixels';
        screensize = num2cell(get(groot, 'MonitorPositions'));

        if monitor > size(screensize, CONST.COLUMNDIR)
            monitor = 1;
        end
        [x, y, w, h] = screensize{monitor, :};

        borderHeight = 84;
        taskBarHeight = 40;
        screenHeight = h - taskBarHeight;

        windowWidth = w/numCols;
        windowHeight = screenHeight/numRows;

        F.Position = [x + windowWidth*(colIndex-1), ...
            y + h - windowHeight*(rowIndex+1), ...
            windowWidth, ...
            max(windowHeight-borderHeight, 100)];
    end
end

