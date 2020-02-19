function [val] = nk(sys, t)
%Changes the clearance parameter depending on the trial time

if t <= sys.GC.t_nk_change
    val = sys.GC.nK/2;
else
    val = sys.GC.nK/2;
end

end

