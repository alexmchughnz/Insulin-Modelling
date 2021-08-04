function PrintTimeTaken(label, patientSet, sinceTime)
    PrintLine();

    tEnd = toc(sinceTime);
    timeTotal = duration(seconds(tEnd));
    totalString = datestr(timeTotal, 'HH:MM:SS');
    
    timePatient = timeTotal / length(patientSet);
    patientString = datestr(timePatient, 'HH:MM:SS');
    
    fprintf("Time taken for %s: %s (Mean of %s per patient)\n", ...
        label, totalString, patientString);    
    
    PrintLine();    
end

