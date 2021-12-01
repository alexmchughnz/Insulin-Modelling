function [P, insertedAt] = InsertData(P, field, tD, vD)
% Inserts, sorted by time, the data points in arrays t and v to the data
% field specified.

    assert(numel(tD) == numel(vD), "Number of time and value points added must match.")

    % Read.
    time = P.data.(field).time;
    value = P.data.(field).value;
        
    % Modify.
    [time, order] = sort([tD; time]);    
    
    value = [vD; value];
    value = value(order);
    
    insertedAt = order(1:numel(tD));  % Indices at start of 'order' are where new points went.

    % Write.  
    P.data.(field).time = time;
    P.data.(field).value = value;

end

