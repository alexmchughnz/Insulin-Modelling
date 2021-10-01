function P = ScalePatientField(P, fieldName, scale)
    
    path = split(fieldName, '.');

    % Edit value of field.
    fieldValue = getfield(P, path{:});
    fieldValue = fieldValue .* scale;
    P = setfield(P, path{:}, fieldValue);

    % Print message.
    if numel(scale) == 1
        scaleString = sprintf("%.4g", scale);
    else
        scaleString = sprintf("%.2g to %.2g", min(scale), max(scale));
    end
    
    message = sprintf("%s scaled by %s.", fieldName, scaleString);
    PrintStatusUpdate(P, message);
    
end