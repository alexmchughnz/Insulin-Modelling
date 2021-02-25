function SaveResults(P, title, saveStruct)
global CONFIG

FILEFORMAT = "%s %s";

filename = sprintf(FILEFORMAT, P.patientCode, title);
path = fullfile(CONFIG.RESULTPATH, filename);

save(path, '-struct', 'saveStruct')

end

