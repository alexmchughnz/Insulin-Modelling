function BolusFunc = MakeBolusFunction(values, times, periods)
% Returns a function which gives delivery rate as a function of time.
% INPUTS:
%   values - array of delivery amounts
%   times - array of delivery times relative to trial start [min]
%   values - array of delivery periods [min]
% OUTPUT:
%   BolusFunc - function handle @(tArray) giving delivery rates for each
%   time in tArray

deliveryRates = values./periods;

% Gives delivery rate at time == t.
GetDeliveryRate = @(t) dot(deliveryRates, (times<=t & t<(times+periods)));

% Applies DeliveryFunc to each time in tArray.
BolusFunc = @(tArray) arrayfun(GetDeliveryRate, tArray);    
end
