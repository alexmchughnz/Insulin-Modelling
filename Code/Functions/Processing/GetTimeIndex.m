function [n] = GetTimeIndex(tSearch, tArray)
% For each t in tSearch, finds the index n of t in tArray.
% Each n is such that tArray(n) < t < tArray(n+1).

n = zeros(size(tSearch));

t0 = tArray(1);
dt = tArray(2) - t0;
for ii = 1:length(tSearch)
    t = tSearch(ii);    
    n(ii) = round((t-t0)/dt) + 1;    
    
    % Fix for rounding up / down too much.
    n = ConstrainArray(n, min(tArray(:)), max(tArray(:)));
end


end

