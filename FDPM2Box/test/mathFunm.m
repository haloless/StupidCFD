function [b] = mathFunm(a,fun)
%mathFunm: Tensor function of eigenvalues
% Decompose A by
% a*v = v*d
% Compute B by
% b = fun(a) = v * fun(d) * inv(v)
%

% TODO check A symmetric

% A*V = V*D
[v,d] = eig(a);

% extract eigenvalues and compute
s = fun(diag(d));

% not work for equal eigenvalues
% b = v * diag(s) * v';

% work around
b = v * diag(s) * inv(v);




return
end