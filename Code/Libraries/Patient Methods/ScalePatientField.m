function P = ScalePatientField(P, scale, varargin)
    path = varargin;

    % Edit value of field.
    fieldValue = getfield(P, path{:});
    fieldValue = fieldValue .* scale;
    P = setfield(P, path{:}, fieldValue);
    
    % Tag patient code with its modifier.
    varName = path{end};
    
    if numel(scale) == 1
        scaleTag = sprintf("%.2g", scale);
    else
        scaleTag = sprintf("[%.2f to %.2f]", min(scale), max(scale));
    end
    tag = sprintf("(%s x %s)", varName, scaleTag);
    
    P.patientCode = strjoin([P.patientCode, tag]);
    
    message = sprintf("%s scaled by %s", varName, scaleTag);
    PrintStatusUpdate(P, message);
end