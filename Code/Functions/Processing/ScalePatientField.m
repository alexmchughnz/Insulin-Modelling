function P = ScalePatientField(P, scale, varargin)
    path = varargin;

    fieldValue = getfield(P, path{:});
    fieldValue = fieldValue .* scale;
    P = setfield(P, path{:}, fieldValue);
    
    P.patientCode = P.patientCode + sprintf(" (%s x %.2g)", ...
        path{end}, scale);
    
    message = sprintf("%s scaled by %.2f.", path{end}, scale);
    PrintStatusUpdate(P, message);
end