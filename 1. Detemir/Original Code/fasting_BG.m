function [BG] = fasting_BG(sys, t)
%returns the correct fasting glucose level for that day

if t < 1000
    BG = sys.GC.fasting_bg1;
else
    BG = sys.GC.fasting_bg2;
end

end

