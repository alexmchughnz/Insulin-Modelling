function[data_datetime, data_cpep, sim_start, sim_end] = get_cpeptide_data(num, txt)

cpep_col = 6;
date_col = 3;

data_cpep = num(:, cpep_col);
data_datetime = txt(:, date_col);
[data_datetime, data_cpep] = data_reformat(data_datetime, data_cpep);
sim_start = data_datetime(1);
sim_end = data_datetime(end);

end