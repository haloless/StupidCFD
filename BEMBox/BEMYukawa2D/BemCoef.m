%%
%% TODO implement a faster numerical quadrature routine
%% for singular (and near-singular) integrals to replace
%% the MATLAB version
%%

function [ aa,bb ] = BemCoef(xin,yin,jelem, singular)

[aa,bb] = BemYukawaCoef2D(xin,yin,jelem, singular);
% [aa,bb] = BemPotentialCoef2D(xin,yin,jelem);

return
end

function [aa,bb] = BemYukawaCoef2D(xin,yin,jelem,singular)

BemMeshGlobals;

% element
al = dlen(jelem);
j1 = node(1,jelem);
j2 = node(2,jelem);

x1 = y(1,j1);
y1 = y(2,j1);
x2 = y(1,j2);
y2 = y(2,j2);

afun = @(s) YukawaG(xin,yin,(1-s)*x1+s*x2,(1-s)*y1+s*y2,kappa);
if singular
    aa = integral(afun, 0.0,1.0);
else
    aa = GaussInteg(afun, 0.0,1.0);
end
aa = aa * al / (pi*2);

bfun = @(s) YukawaF(xin,yin,(1-s)*x1+s*x2,(1-s)*y1+s*y2, dnorm(1,jelem),dnorm(2,jelem), kappa);
if singular
    bb = integral(bfun, 0.0,1.0);
else
    bb = GaussInteg(bfun, 0.0,1.0);
end
bb = bb * al / (pi*2);


return
end

function [g] = YukawaG(x0,y0,x1,y1,k)
r = sqrt((x1-x0).^2 + (y1-y0).^2);
g = besselk(0,r.*k);
return
end

function [f] = YukawaF(x0,y0,x1,y1,nx,ny,k)
% r = sqrt((x1-x0)^2 + (y1-y0)^2);
small = 1.0e-4;
gp = YukawaG(x0,y0,x1+nx*small,y1+ny*small,k);
gn = YukawaG(x0,y0,x1-nx*small,y1-ny*small,k);
f = (gp-gn)./(small*2);
return
end

function [uint] = GaussInteg(ufunc, s0,s1)
ng = 7;
gp = [0.0000000000000000,0.4058451513773972,-0.4058451513773972,0.7415311855993945,-0.7415311855993945,0.9491079123427585,-0.9491079123427585];
gw = [0.4179591836734694, 0.3818300505051189,0.3818300505051189,0.2797053914892766,0.2797053914892766,0.1294849661688697,0.1294849661688697];

sp = (s1+s0)/2 + (s1-s0)/2 .* gp;
up = ufunc(sp);
uint = sum(up.*gw);
uint = (s1-s0)/2 * uint;

return
end




