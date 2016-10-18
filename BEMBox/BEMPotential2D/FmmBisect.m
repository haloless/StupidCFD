
function [ ielem,nsep ] = FmmBisect(x, ielem,ioff, n, xsep, icomp)

nsep = 1;
if n <= 0; return; end

for ifr = 1:n
    % this is the true index in IELEM array
    idx = ifr + ioff - 1;
    
    if x(icomp,ielem(idx)) <= xsep
        if ifr ~= nsep
            isep = nsep + ioff - 1;
            itmp = ielem(isep);
            ielem(isep) = ielem(idx);
            ielem(idx) = itmp;
        end
        nsep = nsep + 1;
    end
end

return
end


