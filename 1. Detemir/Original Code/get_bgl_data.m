function[data_datetime, data_bgl] = get_bgl_data(num, txt)

bgl_col = 4;
date_col = 3;

data_bgl = num(:, bgl_col);
data_datetime = txt(:, date_col);
[data_datetime, data_bgl] = data_reformat(data_datetime, data_bgl);

end