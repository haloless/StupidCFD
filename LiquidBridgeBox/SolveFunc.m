function [x,conv] = SolveFunc(func,x0,tol)
%% A simple implementation of Newton's method
%%

% initial
x = x0;
f = func(x);
conv = 0;

% check initial converge
r0 = abs(f);
disp(['|r0|=',num2str(r0)]);
if r0 <= tol
	conv = 1;
	return;
end

% finite difference to estimate dF/dx
% this step size dh=1.e-6 is roughly (epsilon~1e-16)^(1/3)
dh = 1.0e-6;

%
maxiter = 25;
for iter = 1:maxiter
	df = EstimDeriv(func,x,dh);
	if df == 0
		x = x + dx;
		continue;
	else
		xincr = -f/df;
		x = x + xincr;
	end
	
	f = func(x);
	r = abs(f);
	disp(['iter=',int2str(iter),';|resid|=',num2str(f)]);
	
	% check convergence
	if r <= tol
		conv = 1;
		break;
	end
end

if conv
	disp(['Converged']);
else
	error('Failed');
end



return
end


function [df] = EstimDeriv(func,x,dx)
	fp = func(x+dx);
	fm = func(x-dx);
	df = (fp-fm)/(dx*2);
	return
end


