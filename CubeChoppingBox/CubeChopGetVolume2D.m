
function [ vol ] = CubeChopGetVolume2D(nx,ny,d, lx,ly)
% (x-xc).n = d

xc = lx / 2;
yc = ly / 2;

% flip the orientation
if nx >= 0
    m1 = nx;
else
    m1 = -nx;
end
if ny >= 0
    m2 = ny;
else
    m2 = -ny;
end

alpha = d + xc*m1 + yc*m2;

m1 = m1 * lx;
m2 = m2 * ly;

% put into standard form
msum = m1 + m2;
m1 = m1 / msum;
m2 = m2 / msum;
alpha = alpha / msum;

if m1 > m2
    tmp = m1;
    m1 = m2;
    m2 = tmp;
end

vol = StdChopGetVolume2D(m1,m2,alpha);



return
end

