function P = AddTrialInputs(P)
% Adds bolus arrays to a patient struct.
% Bolus array gives delivery rate over time.

% Insulin Bolus
if isfield(P.data, "vIBolus")
    P.results.IBolus = MakeArray(P, P.data.vIBolus, P.data.tIBolus, P.data.TIBolus);
end

% Glucose Bolus
if isfield(P.data, "vGBolus")
    P.results.GBolus = MakeArray(P, P.data.vGBolus, P.data.tGBolus, P.data.TGBolus);
end

% Exogenous Insulin Function
P.results.Uex = @GetExogenousInsulinInput;
[~, P] = P.results.Uex(P, true);  % To simulate QLocal. Now doesn't need to be done again unless parameters change!

message = "Trial inputs updated.";
PrintStatusUpdate(P, message);

end


function bolusArray = MakeArray(P, values, times, periods)
tArray = P.results.tArray;

deliveryRates = values./periods;

% Gives delivery rate at time == t.
GetDeliveryRate = @(t) dot(deliveryRates, (times<=t & t<(times+periods)));

% Applies DeliveryFunc to each time in tArray.
bolusArray = arrayfun(GetDeliveryRate, tArray);
end