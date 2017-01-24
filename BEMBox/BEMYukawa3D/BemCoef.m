%%
%% TODO implement a faster numerical quadrature routine
%% for singular (and near-singular) integrals to replace
%% the MATLAB version
%%

function [ aa,bb ] = BemCoef(ipos,jelem, singular)

[aa,bb] = BemYukawaCoef3D(ipos,jelem, singular);

return
end



function [aa,bb] = BemYukawaCoef3D(ipos,jelem,singular)

BemMeshGlobals;

% element
ds = facearea(jelem);

j1 = face(jelem,1);
j2 = face(jelem,2);
j3 = face(jelem,3);

x1 = vert(1,j1);
y1 = vert(2,j1);
z1 = vert(3,j1);
x2 = vert(1,j2);
y2 = vert(2,j2);
z2 = vert(3,j2);
x3 = vert(1,j3);
y3 = vert(2,j3);
z3 = vert(3,j3);

nvx = facenvec(1,jelem);
nvy = facenvec(2,jelem);
nvz = facenvec(3,jelem);

xx = ipos(1);
yy = ipos(2);
zz = ipos(3);

afun = @(s,t) YukawaG(xx,yy,zz, ...
(1-s-t)*x1+s*x2+t*x3,...
(1-s-t)*y1+s*y2+t*y3,...
(1-s-t)*z1+s*z2+t*z3,...
kappa);
tmax = @(s) 1-s;

aa = integral2(afun, 0,1, 0,tmax);
aa = aa * ds/0.5;


% if singular
    % aa = integral(afun, 0.0,1.0);
% else
    % aa = GaussInteg(afun, 0.0,1.0);
% end
% aa = aa * al / (pi*2);

bfun = @(s,t) YukawaF(xx,yy,zz,...
(1-s-t)*x1+s*x2+t*x3,...
(1-s-t)*y1+s*y2+t*y3,...
(1-s-t)*z1+s*z2+t*z3,... 
nvx,nvy,nvz,...
kappa);
tmax = @(s) 1-s;

bb = integral2(bfun, 0,1, 0,tmax);
bb = bb * ds/0.5;

% if singular
    % bb = integral(bfun, 0.0,1.0);
% else
    % bb = GaussInteg(bfun, 0.0,1.0);
% end
% bb = bb * al / (pi*2);


return
end

function [g] = YukawaG(x0,y0,z0,x1,y1,z1,k)
r = sqrt((x1-x0).^2 + (y1-y0).^2 + (z1-z0).^2);
g = -exp(-k*r) ./ (4*pi*r);
return
end

function [f] = YukawaF(x0,y0,z0,x1,y1,z1,nx,ny,nz,k)
small = 1.0e-4;
gp = YukawaG(x0,y0,z0,x1+nx*small,y1+ny*small,z1+nz*small,k);
gn = YukawaG(x0,y0,z0,x1-nx*small,y1-ny*small,z1-nz*small,k);
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




