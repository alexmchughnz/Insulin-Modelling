function [P, F] = AddFigure(P, figTitle)
% Adds a figure to patient.
% INPUTS:
%   P        - patient struct
%   figTitle - title of figure
% OUTPUT:
%   P  - patient struct
%   F  - figure handle for debug plot


stack = dbstack(1);
originFunction = stack(1).name;


% Make figure.
fignum = double(string(P.patientNum) + "0" + string(numel(P.figures)));

F = figure(fignum);
F.Name = sprintf("%s: %s", P.patientCode, figTitle);
F.Tag = originFunction;
% F.NumberTitle = "off";

P.figures{end+1} = F;

axes
hold on

end
