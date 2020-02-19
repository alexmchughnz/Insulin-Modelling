function [data_datetime, data] = data_reformat(raw_datetime, raw_data)
%returns the reformatted column - gets rid of the NaN and zero elements
%returns a datetime array where the indexes match the elements in the data
%colum

%setting NaN in cgm data to zero
raw_data(isnan(raw_data)) = 0;

%retrieving indices of all nonzero elements
indices = find(raw_data);

%retrieving datetimes of all nonzero cgm elements
raw_datetime = raw_datetime(indices);

%removing all nonzero elements from cgm column
raw_data(raw_data == 0) = [];
data = raw_data;

%reformating datetime
ts = strrep(raw_datetime,'a.m.','AM');
ts = strrep(ts,'p.m.','PM');
for ii = 1:length(ts)
    try
    % (ADM 19/02/20)
    % Removed 'AM specifier, not present in my data. Catch statement was
    % stripping time data.
    temp = datenum(ts(ii),'dd/mm/yyyy HH:MM:SS'); 
    catch
        temp = datenum(ts(ii),'dd/mm/yyyy');
    end
    data_datetime(ii,1) = datetime(datevec(temp));
end
% data_datetime = datetime(raw_datetime, 'InputFormat', 'dd/MM/yyyy hh:mm:ss aa');

%debugging the .xls midnight NaT error - assumes that there are no two NaT
%are sequential
if length(data_datetime) > 1
    last_date = data_datetime(end-1);
    last_date_string = datestr(last_date);
    just_the_date = split(last_date_string);
    just_the_date = just_the_date{1};
    midnight = just_the_date + " 00:00:00 AM";
    midnight_as_datetime = datetime(midnight, 'InputFormat', 'dd-MMM-yyyy hh:mm:ss aa');
    
    data_datetime(isnat(data_datetime)) = midnight_as_datetime;
end

end

