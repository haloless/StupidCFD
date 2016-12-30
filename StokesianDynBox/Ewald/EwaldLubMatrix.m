function [mat] = EwaldLubMatrix(ewald,np,xp)

lx = ewald.lx;
ly = ewald.ly;
lz = ewald.lz;

% currently FT version
ver = 6;

mat = zeros(np*ver,np*ver);

lubmax = 8.0;

inum = ceil(lubmax/lx);
jnum = ceil(lubmax/ly);
knum = ceil(lubmax/lz);

for i = 1:np
for j = i:np
    % periodic, nearest images
    for ix = -inum:inum
    for iy = -jnum:jnum
    for iz = -knum:knum
        ipos = xp(:,i);
        jpos = xp(:,j) + [lx*ix;ly*iy;lz*iz];
        
        if check_lub(ipos,jpos,lubmax) == 1
            xx = jpos(1) - ipos(1);
            yy = jpos(2) - ipos(2);
            zz = jpos(3) - ipos(3);
            
            r2 = xx^2 + yy^2 + zz^2;
            rr = sqrt(r2);
            
            ex = xx / rr;
            ey = yy / rr;
            ez = zz / rr;
            
            res2b = TwoBodyScalarRes(rr);
            resinf = MinvScalarRes(rr);
            reslub = res2b - resinf;
            % reslub = res2b;
            
            xa11 = reslub(1);
            xa12 = reslub(2);
            ya11 = reslub(3);
            ya12 = reslub(4);
            yb11 = reslub(5);
            yb12 = reslub(6);
            xc11 = reslub(7);
            xc12 = reslub(8);
            yc11 = reslub(9);
            yc12 = reslub(10);
            xg11 = reslub(11);
            xg12 = reslub(12);
            yg11 = reslub(13);
            yg12 = reslub(14);
            yh11 = reslub(15);
            yh12 = reslub(16);
            xm11 = reslub(17);
            xm12 = reslub(18);
            ym11 = reslub(19);
            ym12 = reslub(20);
            zm11 = reslub(21);
            zm12 = reslub(22);
            
            mat = AddMobMatrix(mat,np,i,i,ex,ey,ez, ...
            xa11,ya11, yb11, xc11,yc11, xg11,yg11, yh11, xm11,ym11,zm11);
            
            mat = AddMobMatrix(mat,np,i,j,ex,ey,ez, ...
            xa12,ya12, yb12, xc12,yc12, xg12,yg12, yh12, xm12,ym12,zm12);
            
            mat = AddMobMatrix(mat,np,j,j,-ex,-ey,-ez, ...
            xa11,ya11, yb11, xc11,yc11, xg11,yg11, yh11, xm11,ym11,zm11);
            
            mat = AddMobMatrix(mat,np,j,i,-ex,-ey,-ez, ...
            xa12,ya12, yb12, xc12,yc12, xg12,yg12, yh12, xm12,ym12,zm12);
        end
    end
    end
    end
end
end




return
end

function [ok] = check_lub(x1,x2,lubmax)
r = norm(x1-x2);

if r == 0
    ok = 0;
elseif r >= lubmax
    ok = 0;
else
    ok = 1;
end

return
end


