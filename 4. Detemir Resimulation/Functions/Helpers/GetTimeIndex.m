function [n] = GetTimeIndex(tSearch, tArray)
% For each t in tSearch, finds the index n of t in tArray.
% Each n is such that tArray(n) < t < tArray(n+1).

n = zeros(size(tSearch));

t0 = tArray(1);
dt = tArray(2) - t0;
for ii = 1:length(tSearch)
    t = tSearch(ii);    
    n(ii) = floor((t-t0)/dt) + 1;
end

end

