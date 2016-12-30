% body-centered

clear all;

a = 1.0;

np = 2;

% frac = 0.000125
% frac = 0.008
% frac = 0.027
% frac = 0.064
frac = 0.125

len = (4/3*pi*a^3 * np / frac)^(1/3)

tol = 1.0e-8;
xi = 2 / len;

ewald = EwaldInit(tol, xi, len,len,len)

xp = zeros(3,np);
xp(:,1) = [0,0,0]';
xp(:,2) = [-len/2,-len/2,-len/2]';

up = zeros(6,np);
% particle move U=1
up(1,:) = 1.0;

mob = EwaldMobMatrix(ewald,np,xp);
lub = EwaldLubMatrix(ewald,np,xp);

res = inv(mob);
res = res + lub;

uvec = reshape(up,[],1);

% F = RU
fvec = res * uvec;

fvec(1)










