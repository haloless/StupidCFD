
clear;

k = 1;
% eta = 0.1;
eta = 0.0;
delta = 1.0e-2;

% smoothsign = @(y) y ./ sqrt(y.^2+delta^2));
smoothsign = @(y) tanh(y./delta);

afun = @(t,y) -k*smoothsign(y) - eta*y;

y0 = 1.0;
tmax = 2.0;
% y0 = -1.0e-3;
% tmax = 0.1;
[t,y] = ode45(afun, [0,tmax], [y0]);


%
wfun = @(w) w + delta.*log(tanh(w./delta) ./ (tanh(w./delta)+1));
ts = linspace(0,tmax,11);
yana = [];
w0 = wfun(y0);
for tt = ts
    yana(end+1) = fsolve(@(w) (wfun(w) + k*tt - w0), y0);
end
yana = real(yana);

figure;
plot(t,y,'-x', ts,yana,'+');




