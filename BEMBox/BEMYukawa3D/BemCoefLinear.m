%%
%% TODO implement a faster numerical quadrature routine
%% for singular (and near-singular) integrals to replace
%% the MATLAB version
%%

function [ acoef,bcoef ] = BemCoefLinear(ipos,jelem,v1,v2,v3, level)

BemMeshGlobals;

% element
ds = facearea(jelem);

x1 = v1(1);
y1 = v1(2);
z1 = v1(3);
x2 = v2(1);
y2 = v2(2);
z2 = v2(3);
x3 = v3(1);
y3 = v3(2);
z3 = v3(3);

nvx = facenvec(1,jelem);
nvy = facenvec(2,jelem);
nvz = facenvec(3,jelem);

% vcen = facecent(:,jelem);

% target position
xx = ipos(1);
yy = ipos(2);
zz = ipos(3);

acoef = zeros(3,1);
bcoef = zeros(3,1);

if (level == 1)
	% gauss on triangle
	ng = 7;
	[wg,xg,yg] = TriangleGaussRule(ng);
	
	%
	afun = @(x,y,z) YukawaG(xx,yy,zz, x,y,z, kappa);
	
	aa = LinInteg(afun, v1,v2,v3, ng,wg,xg,yg);
    acoef = acoef + aa;
	
	%
	bfun = @(x,y,z) YukawaF(xx,yy,zz, x,y,z, nvx,nvy,nvz, kappa);
	
	bb = LinInteg(bfun, v1,v2,v3, ng,wg,xg,yg);
    bcoef = bcoef + bb;
	
elseif (level==3)
	% radial singular
	
	ng = 7;
	[wg,xg] = LineGaussRule(ng);
	
	%
	afun = @(x,y,z) YukawaG(xx,yy,zz, x,y,z, kappa);
	
	aa = LinRadInteg(afun, v1,v2,v3, ng,wg,xg,ng,wg,xg);
	acoef = acoef + aa;
	
	%
	bfun = @(x,y,z) YukawaF(xx,yy,zz, x,y,z, nvx,nvy,nvz, kappa);
	
	bb = LinRadInteg(bfun, v1,v2,v3, ng,wg,xg,ng,wg,xg);
	bcoef = bcoef + bb;
	
elseif (level==99)
	% use matlab INTEGRAL2, slow
	smax = @(t) 1-t;
	tmax = @(s) 1-s;
	
	afun = @(s,t) YukawaG(xx,yy,zz, ...
		(1-s-t)*x1+s*x2+t*x3,...
		(1-s-t)*y1+s*y2+t*y3,...
		(1-s-t)*z1+s*z2+t*z3,...
		kappa) .* (1-s-t);
	acoef(1) = integral2(afun, 0,1, 0,tmax);
	
	afun = @(s,t) YukawaG(xx,yy,zz, ...
		(1-s-t)*x1+s*x2+t*x3,...
		(1-s-t)*y1+s*y2+t*y3,...
		(1-s-t)*z1+s*z2+t*z3,...
		kappa) .* s;
	acoef(2) = integral2(afun, 0,1, 0,tmax);
	
	afun = @(s,t) YukawaG(xx,yy,zz, ...
		(1-s-t)*x1+s*x2+t*x3,...
		(1-s-t)*y1+s*y2+t*y3,...
		(1-s-t)*z1+s*z2+t*z3,...
		kappa) .* t;
	acoef(3) = integral2(afun, 0,1, 0,tmax);
	
	
	bfun = @(s,t) YukawaF(xx,yy,zz, ...
		(1-s-t)*x1+s*x2+t*x3,...
		(1-s-t)*y1+s*y2+t*y3,...
		(1-s-t)*z1+s*z2+t*z3,... 
		nvx,nvy,nvz,...
		kappa);
	bcoef(1) = integral2(bfun, 0,1, 0,tmax);
	
	bfun = @(s,t) YukawaF(xx,yy,zz, ...
		(1-s-t)*x1+s*x2+t*x3,...
		(1-s-t)*y1+s*y2+t*y3,...
		(1-s-t)*z1+s*z2+t*z3,... 
		nvx,nvy,nvz,...
		kappa) .* (s);
	bcoef(2) = integral2(bfun, 0,1, 0,tmax);
	
	bfun = @(s,t) YukawaF(xx,yy,zz, ...
		(1-s-t)*x1+s*x2+t*x3,...
		(1-s-t)*y1+s*y2+t*y3,...
		(1-s-t)*z1+s*z2+t*z3,... 
		nvx,nvy,nvz,...
		kappa) .* (t);
	bcoef(3) = integral2(bfun, 0,1, 0,tmax);
	
	acoef = acoef .* ds/0.5;
	bcoef = bcoef .* ds/0.5;
else
	error(['Unknown quadrature level=',int2str(level)]);
end





return
end


function [g] = YukawaG(x0,y0,z0,x1,y1,z1,k)
r = sqrt((x1-x0).^2 + (y1-y0).^2 + (z1-z0).^2);
g = exp(-k*r) ./ (4*pi*r);
g = -g;
return
end

function [f] = YukawaF(x0,y0,z0,x1,y1,z1,nx,ny,nz,k)
% small = 1.0e-4;
% gp = YukawaG(x0,y0,z0,x1+nx*small,y1+ny*small,z1+nz*small,k);
% gn = YukawaG(x0,y0,z0,x1-nx*small,y1-ny*small,z1-nz*small,k);
% f = (gp-gn)./(small*2);
rx = x1 - x0;
ry = y1 - y0;
rz = z1 - z0;
r = sqrt(rx.^2 + ry.^2 + rz.^2) + 1.0e-6;
rinv = 1.0 ./ r;

g = YukawaG(x0,y0,z0,x1,y1,z1,k);
gr = -(rinv+k) .* g;

nr = rx.*nx + ry.*ny + rz.*nz;
f = gr .* rinv .* nr;
return
end


function [cint] = LinInteg(qfunc,x0,x1,x2, ng,wg,xg,yg)

cint = zeros(3,1);

zg = 1-xg-yg;
qx = zg.*x0(1) + xg.*x1(1) + yg.*x2(1);
qy = zg.*x0(2) + xg.*x1(2) + yg.*x2(2);
qz = zg.*x0(3) + xg.*x1(3) + yg.*x2(3);
qval = qfunc(qx,qy,qz);
% qint = sum(qval.*wg);

for ig = 1:ng
    tmp = wg(ig) * qval(ig);
    cint(1) = cint(1) + tmp*zg(ig);
    cint(2) = cint(2) + tmp*xg(ig);
    cint(3) = cint(3) + tmp*yg(ig);
end

% triangle area
aint = cross(x1-x0,x2-x0);
aint = norm(aint);
aint = 0.5 * aint;

cint = cint * aint;


return
end

function [cint] = LinRadInteg(qfunc,x0,x1,x2, nga,wga,xga, ngr,wgr,xgr)

h = 1/sqrt(2);

cint = zeros(3,1);
for i = 1:nga
	% angle
	phi = pi/4 * xga(i);
	cphi = cos(phi);
	sphi = sin(phi);
	cang = cos(phi+pi/4);
	sang = sin(phi+pi/4);
	wphi = wga(i);
	
	for j = 1:ngr
		% radius
		xi = (xgr(j)+1)/2;
		rad = xi * h/cphi;
		wrad = wgr(j);
		
		e1 = rad * cang;
		e2 = rad * sang;
		qpos = x0 + e1*(x1-x0) + e2*(x2-x0);
		
		% function value
		qval = qfunc(qpos(1),qpos(2),qpos(3));
		
		jacob = rad * wphi*wrad * (pi/4) * (0.5*h/cphi);
		
		cint(1) = cint(1) + qval*jacob * (1-e1-e2);
		cint(2) = cint(2) + qval*jacob * (e1);
		cint(3) = cint(3) + qval*jacob * (e2);
	end
end

% patch area
% note coef 0.5 is already treated in preceding integration
aint = cross(x1-x0,x2-x0);
% aint = 0.5 * aint(3);
% aint = aint(3);
aint = norm(aint);

cint = cint .* aint;


return
end



