function output = RoundToMultiple(input, factor, direction)
    % Rounds input to the nearest multiple of factor.
    
    if ~exist('direction', 'var')
        direction = 0;
    end
    
    if direction == 0
        roundFunc = @round;
    elseif direction > 0
        roundFunc = @ceil;
    elseif direction < 0
        roundFunc = @floor;
    end
    
    output = roundFunc(input/factor)*factor;
end

