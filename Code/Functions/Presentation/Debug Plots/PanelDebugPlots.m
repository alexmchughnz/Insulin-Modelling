function PanelDebugPlots(monitor)
% Takes all created debug plot objects and displays them in a grid.
% INPUTS:
%   monitor - integer representing which monitor to display plots on

CONST = LoadConstants();

figData = DebugPlots().FIGURES;

if ~exist('monitor', 'var')
    monitor = 1;
end

GetFigNum = @(figrow) dot([1000 100 1], figrow);

if ~isempty(figData)
    patientCounts = figData(:, 1);
    maxFigs = sum(patientCounts == mode(patientCounts));
    numPatients = length(unique(patientCounts));
    
    seenPatients = [];
    nCounts = zeros(1, numPatients);
    
    for ii = 1 : size(figData, 1)  
        % Grab patient info for figure.
        figRow = figData(ii, :);
        currentPatient = figRow(1);    
        if ~ismember(currentPatient, seenPatients)
            seenPatients = horzcat(seenPatients, currentPatient);
        end
        
        % Make figure.
        figNum = GetFigNum(figRow);
        F = figure(figNum);
        
        % Get position indices.
        pIndex = find(seenPatients == currentPatient);
        nIndex = nCounts(pIndex);
        
        % Position figure.   
        F.Units = 'pixels';     
        screensize = num2cell(get(groot, 'MonitorPositions'));
        
        if monitor > size(screensize, CONST.ROWWISE)
            monitor = 1;
        end
        [x, y, w, h] = screensize{monitor, :};
            

        borderHeight = 84;
        taskBarHeight = 40;
        
        screenHeight = h - taskBarHeight;
        
        windowWidth = w/numPatients;
        windowHeight = screenHeight/maxFigs;
        
        F.Position = [x + windowWidth*(pIndex-1), ...
            y + h - windowHeight*(nIndex+1), ...
            windowWidth, ...
            windowHeight-borderHeight];
        
        % Update pointers.
        nCounts(pIndex) = nCounts(pIndex) + 1;
    end
end
end
