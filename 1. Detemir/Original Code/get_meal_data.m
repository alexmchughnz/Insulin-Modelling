function [sys] = get_meal_data(num, txt, sys)
%returns the meal info for a given patient data
%inputs - num and txt matrices from the patient spreadsheet
%outputs - meal start, duration, carbs (g), glucose (g), and more stuff to
%come...
%assumes trial begins at the first blood reading

meals = num(:, 7);
meals(isnan(meals)) = 0;
meal_indices = find(meals);

times = txt(2:end, 3);
ts = strrep(times,'a.m.','AM');
ts = strrep(ts,'p.m.','PM');
for ii = 1:length(ts)
    try
    temp = datenum(ts(ii),'dd/mm/yyyy HH:MM:SS AM');
    catch
        temp = datenum(ts(ii),'dd/mm/yyyy');
    end
    meal_times(ii,1) = datetime(datevec(temp));
end
% meal_times = datetime(times(meal_indices-1));
sys.Data.meal_start = minutes(meal_times - sys.sim_start_t);

meal_info = num( 2:end, 12:13);
meal_info(isnan(meal_info)) = 0;
carbs = meal_info(:, 1);

meals(meals == 0) = [];
if isempty(meals)
    meals = 20*ones(sum(carbs>0),1);%assume 20 min
end
meal_durations = minutes(meals);
sys.Data.meal_durations = minutes(meal_durations);



sugar = meal_info(:, 2);

sys.Data.meal_start (carbs == 0) = [];
carbs(carbs == 0) = [];

sugar(sugar == 0) = [];

sys.Data.carbs = carbs;
sys.Data.sugar = sugar;

end