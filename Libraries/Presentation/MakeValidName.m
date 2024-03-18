function outString = MakeValidName(inString)
    outString = matlab.lang.makeValidName(inString);
                                      
    % Workaround for x appended to front.
    hasXNow = (outString(1) == 'x');
    hadXOriginally = (inString(1) == 'x');
    
    if hasXNow && ~hadXOriginally
        outString = outString(2:end);
    end
    
end