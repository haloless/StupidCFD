
function [ a ] = FmmUpward(u, a)

BemMeshGlobals;
FmmTreeGlobals;


a(1:nexp+1,1:level(lowlev+2)-1) = 0.0;

for lev = lowlev:-1:2
    ndivx = 2^lev;
    dx = (xmax-xmin)/ndivx;
    dy = (ymax-ymin)/ndivx;
    
    for icell = level(lev+1):level(lev+2)-1
        itr = itree(icell);
        itrx = mod(itr,ndivx);
        itry = floor(itr/ndivx);
        cx = xmin + (itrx+0.5)*dx;
        cy = ymin + (itry+0.5)*dy;
        
        % compute moment for leaf
        if leafflag(icell) == 1
            ioff = loct(icell);
            num = numt(icell);
            uoff = loct(icell);
            acell = FmmMoment(y,node, ielem,ioff, num, nexp, cx,cy, u,uoff, bc,dnorm);
            a(:,icell) = a(:,icell) + acell;
        end
        
        % M2M translation to parent cell
        if lev > 2
            itrxp = floor(itrx/2);
            itryp = floor(itry/2);
            cxp = xmin + (itrxp+0.5)*dx*2;
            cyp = ymin + (itryp+0.5)*dy*2;
            z0 = complex(cx-cxp,cy-cyp);
            io = ifath(icell);
            
            zi = 1.0 + 0.0i;
            for k = 0:nexp
                for m = k:nexp
                    a(m+1,io) = a(m+1,io) + zi*a(m-k+1,icell);
                end
                zi = zi*z0/(k+1);
            end
        end
    end
end

return
end
