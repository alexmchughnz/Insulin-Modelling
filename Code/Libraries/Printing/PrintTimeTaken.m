function PrintTimeTaken(Trial, label)
    PrintLine();

    tEnd = toc(Trial.startTime);
    timeTotal = duration(seconds(tEnd));
    totalString = datestr(timeTotal, 'HH:MM:SS');
    
    timePerPatient = timeTotal / length(Trial.patientSet);
    patientString = datestr(timePerPatient, 'HH:MM:SS');
    
    fprintf("Time taken for %s: %s (Mean of %s per patient)\n", ...
        label, totalString, patientString);    
    
    PrintLine();    
end

