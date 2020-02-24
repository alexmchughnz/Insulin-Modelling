function [kempt] = kempt(current_time,y1,y2,sys)
%Calulcates value of kempt when given glucose concentrations in the two
%stomach compartments
%   Calculates rate at which stomach empties into the gut and is a function of
%   the remining glucose in the two stomach compartments relative to the total
%   glucose ingested (D)

number_of_meals = length(sys.Data.meal_start); %works out how many meals there are
meal_end = sys.Data.meal_start + sys.Data.meal_durations; %end of the meal
%(ADM 24/02/20) Changed 'duration' to 'sys.Data.meal_durations'.

sys.GI.D = 1e-3; %set to a small number by default, then it is overridden by the loop below if reqd

for ii = 1 : number_of_meals
    if current_time > sys.Data.meal_start(ii) && current_time < meal_end(ii)
        sys.GI.D = sys.Data.sugar(ii)/180.156*1000;
    end    
end

b = 0.69;
c = 0.17; %points where kempt intersects kmean
alpha = 5/(2*sys.GI.D*(1-b)); 
beta = 5/(2*sys.GI.D*c);
q = y1 + y2; %sum of compartment 1 and 2 glucose

kempt = sys.GI.kmin + (sys.GI.kmax - sys.GI.kmin)/2 * (tanh(alpha*(q-b*sys.GI.D)) - tanh(beta*(q-c*sys.GI.D)) +2);

end

