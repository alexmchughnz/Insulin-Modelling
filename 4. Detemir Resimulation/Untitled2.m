files = dir(ResultsPath("P*_*_*"));

% Loop through each
for id = 1:length(files)
    % Get the file name (minus the extension)
    [~, f] = fileparts(files(id).name);
      % Convert to number
      usidx = find(f=='_', 1);
      f(usidx) = '';
      
          % If numeric, rename
      movefile(ResultsPath(files(id).name), ResultsPath(sprintf('%s.mat', f)));
end