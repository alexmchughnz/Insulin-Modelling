function [datetime_a, data_cgm_a, datetime_b, data_cgm_b] = get_cgm_data(pt_num)
%retrieves and reformats patient cgm data
%input - patient number i.e 1, 2, 3 etc (not 2 at this stage)
%outputs - cgm data and corresponding times

filename_a = pt_num + "a_cgm.xlsx";
filename_b = pt_num + "b_cgm.xlsx";

addpath '\\file\UsersJ$\jkl59\Home\My Documents\408 FYP\patient data';

[num_a, txt_a, ~] = xlsread(filename_a);
[num_b, txt_b, ~] = xlsread(filename_b);

%bgl data is in column 9
%timedate data is in colum 4
%they start on the 13th row

cgm_col = 9;
date_col = 4;

%get bgl info from column 9
data_cgm_a = num_a(13:end, cgm_col);
data_cgm_b = num_b(13:end, cgm_col);

%get date info from column 4
data_datetime_a = txt_a(13:end, date_col);
data_datetime_b = txt_b(13:end, date_col);

[datetime_a, data_cgm_a] = data_reformat(data_datetime_a, data_cgm_a);
[datetime_b, data_cgm_b] = data_reformat(data_datetime_b, data_cgm_b);

end