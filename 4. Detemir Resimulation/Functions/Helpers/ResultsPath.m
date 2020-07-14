function loadname = ResultsPath(filename)

load('config', 'RESULTPATH');

loadname = fullfile(RESULTPATH, filename);

end
