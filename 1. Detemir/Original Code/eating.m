function [input_rate] = eating(sys, t)
%returns the rate of glucose absoprtion for individual meals
%must be passed in as standard numbers - not durations or times (converted
%with get_meal_data.m)

number_of_meals = length(sys.Data.meal_start); %works out how many meals there are
sys.Data.meal_start = sys.Data.meal_start - 10;
meal_end = sys.Data.meal_start + sys.Data.meal_durations; %end of the meal

input_rate = 0; %by default, then it is overridden by the loop below if reqd

%NOTE - this partially takes into account the carbs...
for ii = 1 : number_of_meals
    if t > sys.Data.meal_start(ii) && t < meal_end(ii)
        input_rate = sys.Data.carbs(ii)/sys.Data.meal_durations(ii)/180.156*1000;% + 0.2*sys.Data.carbs(ii)/sys.Data.meal_durations(ii)/340*1000;
    end    
end

end