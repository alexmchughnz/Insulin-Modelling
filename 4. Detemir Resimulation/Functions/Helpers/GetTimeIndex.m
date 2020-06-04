function [n] = GetTimeIndex(t, tArray)
% Finds the index n of t in tArray, such that tArray(n) < t < tArray(n+1).
t0 = tArray(1);
dt = tArray(2) - t0;
n = floor((t-t0)/dt) + 1;
end

