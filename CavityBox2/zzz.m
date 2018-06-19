
clear;

k = -1.0;

f = @(t,u) k .* u;
g = @(t) exp(k*t);


% dt = 0.005;
% dt = 0.01;
% dt = 0.02;
% dt = 0.04;
% dt = 0.05;
% dt = 0.1;
dt = 0.2;

tmax = 1.0;

nstep = round(tmax/dt);


u = 1.0;
t = 0;
data = [t,u,u];
for istep = 1:nstep
    
    u = u + dt * f(t,u);
    t = t + dt;
    
    data(end+1,:) = [ t, u, g(t) ];
end








