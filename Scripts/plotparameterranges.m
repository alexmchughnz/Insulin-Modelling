clc
clear
close all

% Load data.
patientNums = 'best';
source = "DISST";
setA = GetGridData(source, patientNums);
numA = numel(setA);

patientNums = 'best';
source = "OGTTLui";
setB = GetGridData(source, patientNums);
numB = numel(setB);

patientNums = 'best';
source = "CREBRF";
setC = GetGridData(source, patientNums);
numC = numel(setC);

allPatients = horzcat(setA, struct(), setB, struct(), setC);
numRows = numel(allPatients);

% Plot nL and xL ranges.
figure()
for ii = 1:length(allPatients)
    P = allPatients{ii};
    
    if isempty(fieldnames(P))
        % Hacky to leave a gap between sets.
        continue
    end
    
    xLRange = P.results.optimalxLRange;
    nLRange = P.results.optimalnLRange;
    
    subplot(1,2,1)
    xerrorbar(xLRange, ii, numRows)
    hold on
    
    subplot(1,2,2)
    xerrorbar(nLRange, ii, numRows)
    hold on
end

labelsA = repmat("A", 1, numA) + (numA:-1:1);
labelsB = repmat("B", 1, numB) + (numB:-1:1);
labelsC = repmat("C", 1, numC) + (numC:-1:1);

labels = [labelsC "" labelsB "" labelsA];

subplot(1,2,1)
title("First-pass hepatic insulin clearance")
xlabel("$x_L$ [min$^{-1}$]")
yticks(1:numRows-1)
yticklabels(labels)
ylim([0 numRows])
xlim([0.2 1])
% Add physiological region.
xLPhys = [0.5 0.9];
x = xLPhys([1 1 end end]);
y = ylim;
y = y([1 end end 1]);
patch(x, y, 'r', ...
    'FaceColor', '#D95319', ...
    'FaceAlpha', 0.2, ...
    'EdgeColor', 'none')

subplot(1,2,2)
title("Hepatic insulin clearance")
xlabel("$n_L$")
yticks(1:numRows-1)
yticklabels(labels)
ylim([0 numRows])
xlim([0 0.35])
% Add physiological region.
nLPhys = [0.1 0.3];
x = nLPhys([1 1 end end]);
y = ylim;
y = y([1 end end 1]);
patch(x, y, 'r', ...
    'FaceColor', '#D95319', ...
    'FaceAlpha', 0.2, ...
    'EdgeColor', 'none');


function patients = GetGridData(source, patientNums)
patients = LoadData(source, patientNums);

for ii = 1:length(patients)
    newGrid = false;
    patients{ii} = FindOptimalHepaticClearance(patients{ii}, newGrid);
end
close all
end

function xerrorbar(xrange, ii, nMax)
x = mean(xrange);
y = nMax - ii;
yneg = 0;
ypos = 0;
xneg = abs(x - xrange(1));
xpos = abs(x - xrange(2));
errorbar(x,y,yneg,ypos,xneg,xpos,...
    'Color', [0 1 0],...
    'LineWidth', 1)
end