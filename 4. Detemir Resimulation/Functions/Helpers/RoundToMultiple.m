function output = RoundToMultiple(input, factor)
    % Rounds input to the nearest multiple of factor.
    
    output = round(input/factor)*factor;
end

