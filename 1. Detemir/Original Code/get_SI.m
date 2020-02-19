function [SI] = get_SI(sys, t)
%returns the correct value for SI at corresponding timestep

for ii = 1 : length(sys.GC.SI)
    if t >= ii-1 && t <= ii
        SI = sys.GC.SI(ii);
    end
end

end

