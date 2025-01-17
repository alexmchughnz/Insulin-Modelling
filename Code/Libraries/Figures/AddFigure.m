function [P, F] = AddFigure(P, tag, name)
% Adds a figure to patient.
% INPUTS:
%   P        - patient struct
%   figTitle - title of figure
% OUTPUT:
%   P  - patient struct
%   F  - figure handle for debug plot


F = figure();
F.Name = sprintf("%s: %s", P.patientCode, name);
F.Tag = tag;
F.NumberTitle = "off";

P.figures(end+1) = F;

axes
hold on

end
