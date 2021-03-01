function P = ScalePatientField(scale, P, varargin)
    path = varargin;

    fieldValue = getfield(P, path{:});
    fieldValue = fieldValue .* scale;
    P = setfield(P, path{:}, fieldValue);
    
    fieldName = strjoin(string(path), '.');
    codeParts = split(P.patientCode, " ");
    
    P.patientCode = codeParts(1) + sprintf(" (%s x %.2g)", ...
        fieldName, scale);
    
    message = sprintf("%s scaled by %.2f.", fieldName, scale);
    PrintStatusUpdate(P, message);
end