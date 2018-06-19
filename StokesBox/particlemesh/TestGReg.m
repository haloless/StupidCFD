
clear all;

a = 1.0;

alpha = 0 / a;
xi = 2.0 / a;


rs = linspace(0, 8*a, 101);

gx = [];
gy = [];
gz = [];
for r = rs
	gl = FuncGReg(r,[r,0,0]', alpha,xi);
	gx(end+1) = gl(1,1);
	gy(end+1) = gl(2,2);
	gz(end+1) = gl(3,3);
end

figure;
plot(rs,gx,'x-', rs,gy,'o-');
axis([0,4*a,0,10]);
title(['xi=',num2str(xi),';alpha=',num2str(alpha)]);


