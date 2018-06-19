
clear;

dh = 1.0e-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ADChain.Clear();

% ffunc = @(x1,x2) 2 * x1 .* x2 + 3*sin(x1*4);
ffunc = @(x1,x2) x1./x2 + 2*x1.*x2 + x1*4 + 2*(-x2) - x1.*(x2.^2) + x2.^3 + x1.^(x2/2) + 3*sin(x1*x2) -2*cos(x1*3-x2*1);

x1 = ADScalar(2.0);
x2 = ADScalar(3.0);
f = ffunc(x1,x2);
fval = ffunc(x1.val,x2.val);

[f.val, fval]

ADChain.Propagate();

disp('AutoDiff,NumerDiff');
disp([x1.adj, (ffunc(x1.val+dh,x2.val)-ffunc(x1.val-dh,x2.val))/(dh*2)]);
disp([x2.adj, (ffunc(x1.val,x2.val+dh)-ffunc(x1.val,x2.val-dh))/(dh*2)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ADChain.Clear();

x1 = ADScalar(2.0);
x2 = ADScalar(4.0);
y = -2*x1 + (-3)*x2 + 1;
z = x1 * x2;

% ADChain.Propagate();
% y.propagate();
z.propagate();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ADChain.Clear();

x1 = ADScalar(1.0);
x2 = ADScalar(2.0);
y = atan2(x2,x1);

ADChain.Propagate();
% x1,x2


















