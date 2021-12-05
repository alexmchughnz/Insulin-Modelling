function PanelFigures(Trial, monitor)
% Iterate through all patients' plots and displays them in a grid.

if ~exist('monitor', 'var')
    monitor = 3;
end

CONST = Constants();

numCols = numel(Trial.patientSet);

for pp = 1:numCols
    % Iterate to collect the figures we want to show.
    P = Trial.patientSet{pp};
    allFigs = P.figures;
    displayFigs = {};
    for ii = 1:numel(allFigs)
        F = allFigs(ii);
        tag = F.Tag;
        name = split(F.Name);
        name = name{end};

        try
            toShow = (Trial.figureList.(tag).(name) == true);
        catch err
            if err.identifier == "MATLAB:nonExistentField"
                PrintStatusUpdate(P, "Unknown figure name: " + tag+"."+name)
                toShow = false;
            else
                rethrow(err);
            end
        end

        if toShow
            displayFigs{end+1} = F;
        end
    end

    
    % Position each figure.
    numRows = numel(displayFigs);
    for jj = 1:numRows
        F = figure(displayFigs{jj});

        % Get position indices.
        colIndex = pp;
        rowIndex = jj;

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
            y + h - windowHeight*(rowIndex), ...
            windowWidth, ...
            max(windowHeight-borderHeight, 100)];
    end
end

