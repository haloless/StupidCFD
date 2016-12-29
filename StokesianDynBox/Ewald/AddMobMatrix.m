
function [mat] = AddMobMatrix(mat,n,i,j,ex,ey,ez, ...
xa,ya, yb, xc,yc, ...
xg,yg, yh, xm,ym,zm)

cnt = 6;
nn = n * cnt;
ioff = (i-1) * cnt;
joff = (j-1) * cnt;

ma = MobA(ex,ey,ez, xa,ya);


mb = MobB(ex,ey,ez, yb);

if (i == j)
    mbt = MobB(ex,ey,ez, -yb);
else
    mbt = MobB(ex,ey,ez, yb);
end

mc = MobC(ex,ey,ez, xc,yc);

idx = 1:3;

ii = ioff + idx;
jj = joff + idx;
mat(ii,jj) = mat(ii,jj) + ma;
mat(ii+3,jj) = mat(ii+3,jj) + mb;
mat(ii,jj+3) = mat(ii,jj+3) + mbt;
mat(ii+3,jj+3) = mat(ii+3,jj+3) + mc;


return
end

function [mat] = MobA(ex,ey,ez, xa,ya)
    evec = [ex;ey;ez];
    ee = evec * evec';
    mat = ya.*eye(3) + (xa-ya).*ee;
return
end

function [mat] = MobB(ex,ey,ez, yb)
    b1x = yb * ex;
    b1y = yb * ey;
    b1z = yb * ez;
    
    mat = zeros(3,3);
    mat(1,2) = +b1z;
    mat(1,3) = -b1y;
    mat(2,1) = -b1z;
    mat(2,3) = +b1x;
    mat(3,1) = +b1y;
    mat(3,2) = -b1x;
return
end

function [mat] = MobC(ex,ey,ez, xc,yc)
    evec = [ex;ey;ez];
    ee = evec * evec';
    mat = yc.*eye(3) + (xc-yc).*ee;
return
end

