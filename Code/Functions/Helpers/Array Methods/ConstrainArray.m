function nums = ConstrainArray(nums, lowerB, upperB)

capHigh = @(x) min(x, upperB);
capLow = @(x) max(x, lowerB);

nums = arrayfun(capHigh, nums);
nums = arrayfun(capLow, nums); 

end

