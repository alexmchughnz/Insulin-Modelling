function [utotal] = utotal(t, sys)
%Checks whether the insulin delivery time has been reached and returns
%the rate of change of insulin u_total
%   Detailed explanation goes here
% disp('In u total')
% disp(t)
if t >= (sys.SC.T) && t < (sys.SC.T + 5)
    utotal = sys.SC.Ibolus* 1/5; %rate of change of insulin u_total (assume delivered over half minute period)
else
    utotal = 0;
end

end

