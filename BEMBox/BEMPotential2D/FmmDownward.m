
function [ ax,b ] = FmmDownward(u,ax, a,b)

BemMeshGlobals;
FmmTreeGlobals;

b(1:ntylr+1,1:level(lowlev+2)-1) = 0.0;
ax(1:n) = 0.0;

leaf = 0;
indr = 1;
indi = 1;
for lev = 2:lowlev
    ndivx = 2^lev;
    dx = (xmax-xmin)/ndivx;
    dy = (ymax-ymin)/ndivx;
    
    % loop cells on this level
    for icell = level(lev+1):level(lev+2)-1
        % this cell position
        itr = itree(icell);
        itrx = mod(itr,ndivx);
        itry = floor(itr/ndivx);
        % this cell center
        cx = xmin + (itrx+0.5)*dx;
        cy = ymin + (itry+0.5)*dy;
        % father cell position
        itrxp = floor(itrx/2);
        itryp = floor(itry/2);
        
        % L2L translation from parent level
        if lev > 2
            % father cell center
            cxp = xmin + (itrxp+0.5)*dx*2;
            cyp = ymin + (itryp+0.5)*dy*2;
            z0 = complex(cx-cxp,cy-cyp);
            io = ifath(icell);
            zi = 1.0 + 0.0i;
            for k = 0:ntylr
                for m = 0:ntylr-k
                    b(m+1,icell) = b(m+1,icell) + zi*b(k+m+1,io);
                end
                zi = zi*z0/(k+1);
            end
        end
        
        % interaction list
        for jcell = level(lev+1):level(lev+2)-1
            jtr = itree(jcell);
            jtrx = mod(jtr,ndivx);
            jtry = floor(jtr/ndivx);
            jtrxp = floor(jtrx/2);
            jtryp = floor(jtry/2);
            
            % their parents must be neighbors
            if abs(itrxp-jtrxp)<=1 && abs(itryp-jtryp)<=1
                if abs(itrx-jtrx)>1 || abs(itry-jtry)>1
                    % not direct neighbor, interaction list, use M2L
                    ccx = xmin + (jtrx+0.5)*dx;
                    ccy = ymin + (jtry+0.5)*dy;
                    z0 = complex(cx-ccx,cy-ccy);
                    
                    b(1,icell) = b(1,icell) - log(z0)*a(1,jcell);
                    zo = 1.0;
                    for m = 1:nexp+ntylr
                        zo = zo / z0;
                        kmin = max(0,m-nexp);
                        kmax = min(m,ntylr);
                        sgn = (-1)^kmin;
                        for k = kmin:kmax
                            b(k+1,icell) = b(k+1,icell) + sgn*zo*a(m-k+1,jcell);
                            sgn = -sgn;
                        end
                        zo = zo*m;
                    end
                elseif leafflag(icell)==1 | leafflag(jcell)==1
                    % neighbor leaf, use direct summation
                    if icell == jcell
                        leaf = leaf + 1;
                        leaf3 = leaf*3-1;
                    end
                    
                    ibeg = loct(icell);
                    jbeg = loct(jcell);
                    ax = FmmDirect(ielem,ibeg,jbeg,numt(icell),numt(jcell), ax, u,icell,jcell);
                end
            end
        end
        
        % compute A.x for leaf
        % evaluate local expansion at collocation point
        if leafflag(icell) == 1
            fact = 1.0;
            for itylr = 1:ntylr
                fact = fact / itylr;
                b(itylr+1,icell) = b(itylr+1,icell) * fact;
            end
            for in = 1:numt(icell)
                inax = loct(icell) + in-1;
                indx = ielem(inax);
                zp = b(ntylr+1,icell);
                z0 = complex(x(1,indx)-cx,x(2,indx)-cy);
                for itylr = ntylr-1:-1:0
                    zp = zp*z0 + b(itylr+1,icell);
                end
                zp = zp / (pi*2);
                ax(inax) = ax(inax) + real(zp);
            end
        end
    end
end

return
end


function [ax] = FmmDirect(elem,ibeg,jbeg,ni,nj,ax,u,icell,jcell)

BemMeshGlobals;

for j = 1:nj
    % j element
    % jind = jnod(j);
    jj = jbeg + j - 1;
    jind = elem(jj);
    
    for i = 1:ni
        % iind = inod(i);
        ii = ibeg + i - 1;
        iind = elem(ii);
        
        [aa,bb] = BemElemInteg(iind,jind);
        
        if bc(1,jind) == bc_dir
            ax(ii) = ax(ii) - aa*u(jj);
        elseif bc(1,jind) == bc_neu
            ax(ii) = ax(ii) + bb*u(jj);
        end
    end
end


return
end






