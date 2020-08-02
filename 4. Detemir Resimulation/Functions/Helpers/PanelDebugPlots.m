function F = PanelDebugPlots(monitor)
set(groot, 'defaultFigureUnits', 'pixels');

GetFigNum = @(figrow) dot([100 10 1], figrow);

global DEBUGPLOTS
figData = DEBUGPLOTS.FIGURES;

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
        screensize = num2cell(get(groot, 'MonitorPositions'));
        [x, y, w, h] = screensize{monitor, :};
        
        position = num2cell(F.Position);
        [~, ~, ~, h0] = position{1, :};
        
        F.Position = [x + w*(pIndex-1)/3, ...
            (h-1.2*h0) - h*(nIndex/(maxFigs+1)), ...
            w/numPatients, ...
            h0];
        
        % Update pointers.
        nCounts(pIndex) = nCounts(pIndex) + 1;
    end
end
end
