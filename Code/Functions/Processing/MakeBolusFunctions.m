function P = MakeBolusFunctions(P)
% Adds a bolus function to a patient struct.
% Bolus function gives delivery rate as a function of time.

% Insulin Bolus
if isfield(P.data, "vIBolus")
    P.data.IBolus = MakeFunc(P.data.vIBolus, P.data.tIBolus, P.data.TIBolus);
end

% Glucose Bolus
if isfield(P.data, "vGBolus")
    P.data.GBolus = MakeFunc(P.data.vGBolus, P.data.tGBolus, P.data.TGBolus);
end
end


function bolusFunc = MakeFunc(values, times, periods)
deliveryRates = values./periods;

% Gives delivery rate at time == t.
GetDeliveryRate = @(t) dot(deliveryRates, (times<=t & t<(times+periods)));

% Applies DeliveryFunc to each time in tArray.
bolusFunc = @(tArray) arrayfun(GetDeliveryRate, tArray);
end