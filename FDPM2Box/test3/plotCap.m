function [] = plotCap(par, a, varargin)
%plotCap

beta = par.beta;
M = par.M;
pt = par.ptens;

b = beta;

tt = linspace(0,pi/2,91);
tt = tt + pi/2;

xx = (pt-a) + b*a*cos(tt);
yy = M/sqrt(3)*a*sin(tt);

plot(xx,yy,varargin{:});


return
end


