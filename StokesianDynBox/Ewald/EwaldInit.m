
function [ ewald ] = EwaldInit(tol,xi, lx,ly,lz)

ewald = struct();

%
ewald.eps = tol;
ewald.xi = xi;
ewald.lx = lx;
ewald.ly = ly;
ewald.lz = lz;
ewald.pivol = pi / (lx*ly*lz);
%
% diagonal part
%
xi2 = xi * xi;
xia2 = xi2;
xiaspi = xi / sqrt(pi);

ewald.self_a = 1.0 - xiaspi * (6 - 40/3*xia2);
ewald.self_c = 0.75 - xiaspi * xia2 * 10;
ewald.self_m = 0.9 - xiaspi * xia2 * (12 - 30.24*xia2);

%
% lattice for Ewald Sum
%
red_factor = 0.99;
lmin = min([lx,ly,lz]);
lmax = max([lx,ly,lz]);
dr = lmin;
dk = 2*pi/lmax;

% real space
rmax = 10 * sqrt(-log(tol)) / xi;
while 1
    rtest = rmax * red_factor;
    tmp = EwaldReal(rtest,xi) + 6*EwaldReal(rtest+dr,xi);
    if tmp < tol
        rmax = rtest;
    else
        break;
    end
end
ewald.rmax = rmax;
ewald.rnumx = ceil(rmax/lx);
ewald.rnumy = ceil(rmax/ly);
ewald.rnumz = ceil(rmax/lz);

% wave space
kmax = 10 * sqrt(-log(tol)) * 2*xi;
while 1
    ktest = kmax * red_factor;
    tmp = EwaldRecip(ktest,xi) + 6*EwaldRecip(ktest+dk,xi);
    if tmp < tol
        kmax = ktest;
    else
        break;
    end
end
ewald.kmax = kmax;
ewald.knumx = floor(kmax*lx/(pi*2));
ewald.knumy = floor(kmax*ly/(pi*2));
ewald.knumz = floor(kmax*lz/(pi*2));


return
end


function [val] = EwaldReal(r,xi)

xir = xi * r / 1.15;
xir2 = xir * xir;
r2 = r * r;

val = erfc(xir)*(0.75+1.5/r2)/r + ...
exp(-xir2) * xi/sqrt(pi) * (4.5+3*xir2+(3+xir2*(14+xir2*(4+xir2)))/r2);

return
end

function [val] = EwaldRecip(k,xi)

kx = k / xi / 1.15;
kx2 = kx * kx;
k2 = k * k;

val = 6*pi * exp(-kx2/4) * (1+k2/3)/k2 * (1+kx2*(0.25+kx2/8));

return
end







