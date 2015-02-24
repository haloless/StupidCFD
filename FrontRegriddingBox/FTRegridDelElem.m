
function [ neldel ] = FTRegridDelElem()

FTRegridGlobals;

neldel = 0;

ke = ffe;
% this loop takes MAXEL, because elements are being dynamically deleted
for kk = 1:maxel
    %
    [s,n] = felside(ke,icp,pt);
    
    flag = 0;
    if s(1) < amin
        flag = 1;
    end
    
    kenew = ke;
    if (flag == 1)
        kenew = fdelel(ke,n(1));
        
        neldel = neldel + 2;
    end
    
    if (ke==lfe); break; end
    
    if kenew==ke
        ke = elcon(ke);
    else
        ke = kenew;
    end
end


return
end


