
clear;

funj0 = @(x) besselj(0,x);
funj1 = @(x) besselj(1,x);

% xguess = [4, 6, 10, 12, 16, 20];
% lam = [];
% for xg = xguess
    % lam(end+1) = fsolve(funj1, xg);
% end

% lam = [];
% for i = [4, 6:2:1000]
    % lam(end+1) = fzero(funj1, i);
% end
% lam = uniquetol(lam, 1.0e-5);
% lam = unique(lam);

lam = besselzero(1, 1000, 1);
lam = lam.';


if 0
figure;
xs = 0:0.01:20;
js = funj1(xs);
plot(xs,js,'-', lam,zeros(size(lam)),'x');
end

ncut = numel(lam);

W = 2*pi*63/60;
R = 0.06;
nu = 1.0e-6;
% nu = 1.0e-3;

tstar = R^2 / nu / lam(1)^2

t = 0;
% t = 0.3
t = 200

rs = linspace(0, R, 65);
us = [];
for r = rs
    u = 0;
    for l = lam
        ul = funj1(l*r/R) / (l*funj0(l)) * exp(-l^2*nu*t/R^2);
        u = u + ul;
    end
    u = u * (-2*W*R);
    
    if 0
        % spin-down
        u = u;
    else
        % spin-up, impose rotational speed
        u = (W*r - u);
    end
    
    % normalized by ref vel
    % u = u / (W*R);
    
    us(end+1) = u;
end

uu = us';

if (1)
    figure;
    plot(rs,us,'x-');
end








