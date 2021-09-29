function P = ScalePatientField(P, scale, varargin)
    path = varargin;

    % Edit value of field.
    fieldValue = getfield(P, path{:});
    fieldValue = fieldValue .* scale;
    P = setfield(P, path{:}, fieldValue);
    
    % Tag patient code with its modifier.
    varName = path{end};
    
    if numel(scale) == 1
        scaleTag = sprintf("%.4g", scale);
    else
        scaleTag = sprintf("[%.2f to %.2f]", min(scale), max(scale));
    end
    tag = sprintf("%s * %s", varName, scaleTag);
    
    P = TagPatientCode(P, tag);
end