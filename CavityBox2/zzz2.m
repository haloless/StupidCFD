
clear;

g = 1.0;
v = 1.0;
b = 1.0;



% dt = 0.001;
% dt = 0.002;
% dt = 0.005;
% dt = 0.01;
% dt = 0.02;
% dt = 0.05;
dt = 0.1;
% dt = 0.2;
% dt = 0.5;
% dt = 1.0;



tmax = 1.0;
nstep = round(tmax/dt);

method = 3;

u = 0;
t = 0;
for istep = 1:nstep
	switch method 
	case 0
		% euler explicit
		du = g + b * (v - u);
		u = u + dt * du;
	case 1
		% euler implicit
		u = (u + dt*g + dt*b*v) / (1+dt*b);
	case 2
		% 
		uu = (u + dt*b*v) / (1+dt*b);
		u = uu + dt*g;
	case 3
		% 
		uu = u + dt*b*v;
		u = (uu + dt*g) / (1+dt*b);
	case 4
		% 
		uu = u + dt*g + dt*b*v;
		u = (uu + dt*g) / (1+dt*b);
	otherwise
		error('invalid method')
	end
	
	t = t + dt;
end


dt
uana = (g/b+v) * (1-exp(-b*t));
err = abs(u - uana)






