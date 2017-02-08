function [phi,gphi,hphi] = SDFPotential(sdf,x,y)

xmin = sdf.xmin;
ymin = sdf.ymin;
dh = sdf.dh;

% locate position in SDF grid
ibase = floor((x-xmin)/dh) + 1;
jbase = floor((y-ymin)/dh) + 1;
if ibase<1 || ibase>=sdf.nx || jbase<1 || jbase>=sdf.ny
	error(['position overflow']);
end

xbase = (ibase-1)*dh + xmin;
ybase = (jbase-1)*dh + ymin;

xx = (x-xbase) / dh;
yy = (y-ybase) / dh;
if xx<0 || xx>1 || yy<0 || yy>1
	error(['interpolation error']);
end

% local data for interpolation
data = sdf.phig(ibase-1:ibase+2,jbase-1:jbase+2);

% interpolation coef.
[u,ux,uxx] = CubicInterpCoef(xx,dh);
[v,vy,vyy] = CubicInterpCoef(yy,dh);

phi = 0;
phix = 0;
phiy = 0;
phixx = 0;
phixy = 0;
phiyy = 0;

for ii = 1:4
for jj = 1:4
	phi = phi + u(ii)*v(jj)*data(ii,jj);
	
	phix = phix + ux(ii)*v(jj)*data(ii,jj);
	phiy = phiy + u(ii)*vy(jj)*data(ii,jj);
	
	phixx = phixx + uxx(ii)*v(jj)*data(ii,jj);
	phixy = phixy + ux(ii)*vy(jj)*data(ii,jj);
	phiyy = phiyy + u(ii)*vyy(jj)*data(ii,jj);
end
end

gphi = [ phix; phiy ];
hphi = [ phixx, phixy; phixy, phiyy ];

return
end

function [c,cd,cdd] = CubicInterpCoef(x,dh)
	x2 = x * x;
	
	c = [ -x*(x-1)*(x-2)/6; (x+1)*(x-1)*(x-2)/2; -(x+1)*x*(x-2)/2; (x+1)*x*(x-1)/6 ];
	
	cd = [ -(3*x2-6*x+2)/6; (3*x2-4*x-1)/2; -(3*x2-2*x-2)/2; (3*x2-1)/6 ];
	cd = cd ./ dh;
	
	cdd = [ 1-x; 3*x-2; 1-3*x; x ];
	cdd = cdd ./ (dh^2);
	
	return
end



