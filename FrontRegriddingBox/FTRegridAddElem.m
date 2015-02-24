
function [ neladd ] = FTRegridAddElem()

FTRegridGlobals;

neladd = 0;

ke = ffe;
% this loop takes MAXEL, because elements are being dynamically added
for kk = 1:maxel 
    %
    [s,n] = felside(ke,icp,pt);
    
    flag = 0;
    if (s(3)>amax)
        flag = 1;
    end
    if (s(1)>amin)
        aspr = faspr(ke,icp,pt);
        if (aspr > aspmax)
            flag = 1;
        end
    end
    
    if (flag == 1) 
        % add 2 elements, one for current, one for its neighbor
        if (np+1>=maxpt); error('point overflow'); end
        if (nelem+2>=maxel); error('element overflow'); end
        
        %
        faddel(ke,n(3));
        
        neladd = neladd + 2;
    end
    
    if (ke==lfe); break; end
    ke = elcon(ke);
end


return
end


