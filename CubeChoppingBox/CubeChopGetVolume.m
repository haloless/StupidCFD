
function [ vol ] = CubeChopGetVolume(nx,ny,nz,d, lx,ly,lz)
% (x-xc).n = d

xc = lx / 2;
yc = ly / 2;
zc = ly / 2;


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
if nz >= 0
    m3 = nz;
else
    m3 = -nz;
end

alpha = d + xc*m1 + yc*m2 + zc*m3;

m1 = m1 * lx;
m2 = m2 * ly;
m3 = m3 * lz;

% put into standard form
msum = m1 + m2 + m3;
m1 = m1 / msum;
m2 = m2 / msum;
m3 = m3 / msum;
alpha = alpha / msum;

if m1 > m2
    tmp = m1;
    m1 = m2;
    m2 = tmp;
end
if m1 > m3
    tmp = m1;
    m1 = m3;
    m3 = tmp;
end
if m2 > m3
    tmp = m2;
    m2 = m3;
    m3 = tmp;
end

vol = StdChopGetVolume(m1,m2,m3,alpha);



return
end

