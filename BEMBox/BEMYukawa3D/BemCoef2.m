%%
%% TODO implement a faster numerical quadrature routine
%% for singular (and near-singular) integrals to replace
%% the MATLAB version
%%

function [ aa,bb ] = BemCoef2(ipos,jelem, level)

BemMeshGlobals;

% ipos
% facecent(:,jelem)

% element
ds = facearea(jelem);

j1 = face(jelem,1);
j2 = face(jelem,2);
j3 = face(jelem,3);

v1 = vert(:,j1);
v2 = vert(:,j2);
v3 = vert(:,j3);
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

vcen = facecent(:,jelem);

% target position
xx = ipos(1);
yy = ipos(2);
zz = ipos(3);


aa = 0;
bb = 0;

if (level == 0)
	% use midpoint value
	% in most cases, do not use this
	aa = YukawaG(xx,yy,zz, vcen(1),vcen(2),vcen(3), kappa) * ds;
	bb = YukawaF(xx,yy,zz, vcen(1),vcen(3),vcen(3), nvx,nvy,nvz, kappa) * ds;
elseif (level == 1)
	% gauss on triangle
	ng = 7;
	[wg,xg,yg] = TriangleGaussRule(ng);
	
	afun = @(x,y,z) YukawaG(xx,yy,zz, x,y,z, kappa);
	
	aa = 0;
	aa = aa + TriangleGaussInteg(afun,v1,v2,v3, ng,wg,xg,yg);
	
	bfun = @(x,y,z) YukawaF(xx,yy,zz, x,y,z, nvx,nvy,nvz, kappa);
	
	bb = 0;
	bb = bb + TriangleGaussInteg(bfun,v1,v2,v3, ng,wg,xg,yg);
	
elseif (level==2)
	% gauss on three sub-triangle
	
	ng = 7;
	[wg,xg,yg] = TriangleGaussRule(ng);
	
	afun = @(x,y,z) YukawaG(xx,yy,zz, x,y,z, kappa);
	
	aa = 0;
	aa = aa + TriangleGaussInteg(afun,vcen,v1,v2, ng,wg,xg,yg);
	aa = aa + TriangleGaussInteg(afun,vcen,v2,v3, ng,wg,xg,yg);
	aa = aa + TriangleGaussInteg(afun,vcen,v3,v1, ng,wg,xg,yg);
	
	bfun = @(x,y,z) YukawaF(xx,yy,zz, x,y,z, nvx,nvy,nvz, kappa);
	
	bb = 0;
	bb = bb + TriangleGaussInteg(bfun,vcen,v1,v2, ng,wg,xg,yg);
	bb = bb + TriangleGaussInteg(bfun,vcen,v2,v3, ng,wg,xg,yg);
	bb = bb + TriangleGaussInteg(bfun,vcen,v3,v1, ng,wg,xg,yg);
	
elseif (level==3)
	% radial singular
	
	ng = 7;
	[wg,xg] = LineGaussRule(ng);
	
	afun = @(x,y,z) YukawaG(xx,yy,zz, x,y,z, kappa);
	v4 = (v1+v2) / 2;
	v5 = (v2+v3) / 2;
	v6 = (v3+v1) / 2;
	% v7 = (v1+v4) / 2;
	% v8 = (v4+v2) / 2;
	% v9 = (v2+v5) / 2;
	% v10 = (v5+v3) / 2;
	% v11 = (v3+v6) / 2;
	% v12 = (v6+v1) / 2;
	
	aa = 0;
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v1,v2, ng,wg,xg,ng,wg,xg);
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v2,v3, ng,wg,xg,ng,wg,xg);
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v3,v1, ng,wg,xg,ng,wg,xg);
	aa = aa + TriangleRadSingularInteg(afun,vcen,v1,v4, ng,wg,xg,ng,wg,xg);
	aa = aa + TriangleRadSingularInteg(afun,vcen,v4,v2, ng,wg,xg,ng,wg,xg);
	aa = aa + TriangleRadSingularInteg(afun,vcen,v2,v5, ng,wg,xg,ng,wg,xg);
	aa = aa + TriangleRadSingularInteg(afun,vcen,v5,v3, ng,wg,xg,ng,wg,xg);
	aa = aa + TriangleRadSingularInteg(afun,vcen,v3,v6, ng,wg,xg,ng,wg,xg);
	aa = aa + TriangleRadSingularInteg(afun,vcen,v6,v1, ng,wg,xg,ng,wg,xg);
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v1,v7, ng,wg,xg,ng,wg,xg);
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v7,v4, ng,wg,xg,ng,wg,xg);
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v4,v8, ng,wg,xg,ng,wg,xg);
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v8,v2, ng,wg,xg,ng,wg,xg);
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v2,v9, ng,wg,xg,ng,wg,xg);
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v9,v5, ng,wg,xg,ng,wg,xg);
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v5,v10, ng,wg,xg,ng,wg,xg);
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v10,v3, ng,wg,xg,ng,wg,xg);
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v3,v11, ng,wg,xg,ng,wg,xg);
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v11,v6, ng,wg,xg,ng,wg,xg);
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v6,v12, ng,wg,xg,ng,wg,xg);
	% aa = aa + TriangleRadSingularInteg(afun,vcen,v12,v1, ng,wg,xg,ng,wg,xg);
	
	bfun = @(x,y,z) YukawaF(xx,yy,zz, x,y,z, nvx,nvy,nvz, kappa);
	bb = 0;
	bb = bb + TriangleRadSingularInteg(bfun,vcen,v1,v2, ng,wg,xg,ng,wg,xg);
	bb = bb + TriangleRadSingularInteg(bfun,vcen,v2,v3, ng,wg,xg,ng,wg,xg);
	bb = bb + TriangleRadSingularInteg(bfun,vcen,v3,v1, ng,wg,xg,ng,wg,xg);
	
elseif (level==99)
	% use matlab INTEGRAL2, slow
	
	afun = @(s,t) YukawaG(xx,yy,zz, ...
		(1-s-t)*x1+s*x2+t*x3,...
		(1-s-t)*y1+s*y2+t*y3,...
		(1-s-t)*z1+s*z2+t*z3,...
		kappa);
	tmax = @(s) 1-s;

	aa = integral2(afun, 0,1, 0,tmax);
	aa = aa * ds/0.5;

	bfun = @(s,t) YukawaF(xx,yy,zz, ...
		(1-s-t)*x1+s*x2+t*x3,...
		(1-s-t)*y1+s*y2+t*y3,...
		(1-s-t)*z1+s*z2+t*z3,... 
		nvx,nvy,nvz,...
		kappa);
	tmax = @(s) 1-s;

	bb = integral2(bfun, 0,1, 0,tmax);
	bb = bb * ds/0.5;
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
r = sqrt(rx.^2 + ry.^2 + rz.^2);
rinv = 1.0 ./ r;

g = YukawaG(x0,y0,z0,x1,y1,z1,k);
gr = -(rinv+k) .* g;

nr = rx.*nx + ry.*ny + rz.*nz;
f = gr .* rinv .* nr;
return
end

% function [uint] = GaussInteg(ufunc, s0,s1)
% ng = 7;
% gp = [0.0000000000000000,0.4058451513773972,-0.4058451513773972,0.7415311855993945,-0.7415311855993945,0.9491079123427585,-0.9491079123427585];
% gw = [0.4179591836734694, 0.3818300505051189,0.3818300505051189,0.2797053914892766,0.2797053914892766,0.1294849661688697,0.1294849661688697];

% sp = (s1+s0)/2 + (s1-s0)/2 .* gp;
% up = ufunc(sp);
% uint = sum(up.*gw);
% uint = (s1-s0)/2 * uint;

% return
% end




