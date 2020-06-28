function F = PanelDebugPlots()

GetFigNum = @(figrow) dot([100 10 1], figrow);

global DEBUGPLOTS
figData = DEBUGPLOTS.FIGURES;

if ~isempty(figData)
    numPatients = 3;
    patientCounts = figData(:, 1);
    maxFigs = sum(patientCounts == mode(patientCounts));
    
    nCount = 0;
    pCount = 0;
    prevPatient = 0;
    for ii = 1 : size(figData, 1)
        % Get figure.
        figRow = figData(ii, :);
        figNum = GetFigNum(figRow);
        F = figure(figNum);
        
        % Update patient numbers if required.
        currentPatient = figRow(1);
        if currentPatient ~= prevPatient
            pCount = pCount + 1;
            nCount = 0;
        end
        
        % Position figure.
        monitor = 1;
        
        screensize = num2cell(get(groot, 'MonitorPositions'));
        [x, y, w, h] = screensize{monitor, :};
        
        position = num2cell(F.Position);
        [~, ~, ~, h0] = position{1, :};
        
        F.Position = [x + w*(pCount-1)/3, ...
            (h-1.2*h0) - h*(nCount/(maxFigs+1)), ...
            w/numPatients, ...
            h0];
        
        % Update pointers.
        prevPatient = currentPatient;
        nCount = nCount + 1;
    end
end
end
