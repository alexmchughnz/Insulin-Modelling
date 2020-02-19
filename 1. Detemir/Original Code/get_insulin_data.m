function[data_datetime, data_insulin] = get_insulin_data(num, txt)

cpep_col = 5;
date_col = 3;

data_insulin = num(:, cpep_col);
data_datetime = txt(:, date_col);
[data_datetime, data_insulin] = data_reformat(data_datetime, data_insulin);

end