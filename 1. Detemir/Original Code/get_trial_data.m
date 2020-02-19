function [trial_start, trial_end] = get_trial_data(txt)
%gets the first and last trial dates
%hard coded to start at the 4th element in date-time array - must change this if
%the patient data file is different

date_col = 3;
data_datetime = txt(:, date_col);


ts = strrep(data_datetime(2),'a.m.','AM');
ts = strrep(ts,'p.m.','PM');
% (ADM 19/02/20)
% Removed 'AM' specifier, not present in my data. Catch statement was
% stripping time data. Same for te, below.
ts = datenum(ts,'dd/mm/yyyy HH:MM:SS');

te = strrep(data_datetime(end),'a.m.','AM');
te = strrep(te,'p.m.','PM');
te = datenum(te,'dd/mm/yyyy HH:MM:SS');
trial_start = datetime(datevec(ts));
trial_end =datetime(datevec(te));


end
