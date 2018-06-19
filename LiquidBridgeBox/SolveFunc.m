function [x,conv] = SolveFunc(func,x0,tol, varargin)
%% A simple implementation of Newton's method
%%

% finite difference to estimate dF/dx
% this step size dh=1.e-6 is roughly (epsilon~1e-16)^(1/3)
dh = 1.0e-6;

deriv = 0;

kvpairs = varargin;
for i = 1:2:length(kvpairs)
    key = kvpairs{i};
    val = kvpairs{i+1};
    
    switch key
    case 'Derivative'
        switch val
        case 'numer-diff'
            deriv = 0;
        case 'auto-diff'
            deriv = 1;
        otherwise
            error('Unknown Derivative = %s', val);
        end
    otherwise
        error('Unknown %s',kvpairs{i});
    end
end

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


%
maxiter = 25;
for iter = 1:maxiter
    
    if deriv == 0
        df = EstimDeriv(func,x,dh);
    elseif deriv == 1
        df = AutoDeriv(func, x);
    end
    
	if df == 0
		x = x + dh;
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
	% error('Failed');
	warning('Failed');
end



return
end


function [df] = EstimDeriv(func,x,dx)
	fp = func(x+dx);
	fm = func(x-dx);
	df = (fp-fm)/(dx*2);
	return
end

function [df] = AutoDeriv(func,x)
    ADChain.Clear();
    xwrap = ADScalar(x);
    fwrap = func(xwrap);
    ADChain.Propagate();
    df = xwrap.adj;
	return
end


