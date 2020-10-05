patientNums = [2 3 5 7 8 9 13 14 24 25];
source = "DISST";
DISSTPatients = GetData(source, patientNums);

patientNums = [12 128 146 160 166 169 171 196 198 216];
source = "CREBRF";
CREBRFPatients = GetData(source, patientNums);

% Retrieve nL and xL ranges.
figure()


patients = horzcat(DISSTPatients, struct(), CREBRFPatients);
for ii = 1:length(patients)
    P = patients{ii};
    
    if isempty(fieldnames(P))
        % Hacky to leave a gap between A and B.
        continue
    end
    
    xLRange = P.results.optimalxLRange;
    nLRange = P.results.optimalnLRange;
    
    subplot(1,2,1)
    xerrorbar(xLRange, ii)
    hold on
    
    subplot(1,2,2)
    xerrorbar(nLRange, ii)
    hold on
end

labels = [repmat("B", 1, 10) + (10:-1:1), "", repmat("A", 1, 10) + (10:-1:1)];

subplot(1,2,1)
xlabel("$x_L$ [min$^{-1}$]")
yticks(1:21)
yticklabels(labels)
ylim([0 22])
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
xlabel("$n_L$")
yticks(1:21)
yticklabels(labels)
ylim([0 22])
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


function patients = GetData(source, patientNums)
patients = makedata(source, patientNums);

for ii = 1:length(patients)
    patients{ii} = FindOptimalHepaticClearance(patients{ii}, 'load');
end
close all
end

function xerrorbar(xrange, ii)
x = mean(xrange);
y = 22 - ii;
yneg = 0;
ypos = 0;
xneg = abs(x - xrange(1));
xpos = abs(x - xrange(2));
errorbar(x,y,yneg,ypos,xneg,xpos, 'k')
end