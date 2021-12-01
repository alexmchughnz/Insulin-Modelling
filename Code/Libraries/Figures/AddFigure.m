function [P, F] = AddFigure(P, tag, name)
% Adds a figure to patient.
% INPUTS:
%   P        - patient struct
%   figTitle - title of figure
% OUTPUT:
%   P  - patient struct
%   F  - figure handle for debug plot

% Make figure.
fignum = double(string(P.patientNum) + "0" + string(numel(P.figures)));

F = figure(fignum);
F.Name = sprintf("%s: %s", P.patientCode, name);
F.Tag = tag;
F.NumberTitle = "off";

P.figures{end+1} = F;

axes
hold on

end
