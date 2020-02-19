function [datetime_bolus, bolus] = get_bolus_data(num, txt)

bolus_col = 15;
date_col = 3;

bolus = num(:, bolus_col);
datetime_bolus = txt(:, date_col);

[datetime_bolus, bolus] = data_reformat(datetime_bolus, bolus);

end

