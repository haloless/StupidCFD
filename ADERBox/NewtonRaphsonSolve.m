% ## Copyright (C) 2014 homu
% ## 
% ## This program is free software; you can redistribute it and/or modify
% ## it under the terms of the GNU General Public License as published by
% ## the Free Software Foundation; either version 3 of the License, or
% ## (at your option) any later version.
% ## 
% ## This program is distributed in the hope that it will be useful,
% ## but WITHOUT ANY WARRANTY; without even the implied warranty of
% ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% ## GNU General Public License for more details.
% ## 
% ## You should have received a copy of the GNU General Public License
% ## along with Octave; see the file COPYING.  If not, see
% ## <http://www.gnu.org/licenses/>.

% ## NewtonRaphsonSolve

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-01-22

function [ retval, x,resnorm,niter ] = NewtonRaphsonSolve (F, dFdx, x0, options)
% Description

% OPTIONS is a structure of solver options created using OPTIMSET.
% EG: options = optimset('TolX', 0.001).
% The following options can be set:
% * OPTIONS.TOLFUN is the maximum tolerance of the norm of the residuals.
%   [1e-6]
% * OPTIONS.TOLX is the minimum tolerance of the relative maximum stepsize.
%   [1e-6]
% * OPTIONS.MAXITER is the maximum number of iterations before giving up.
%   [100]
% * OPTIONS.DISPLAY sets the level of display: {'off', 'iter'}.
%   ['iter']

if (isrow(x0)); x0 = x0'; end

defopts = optimset('TolX',1e-12, 'TolFun',1e-9, 'MaxIter',100, 'Display','off');
if (nargin < 4)
    options = defopts;
else
    opt1 = defopts;
    options = optimset(opt1, options);
end

tol_x = optimget(options, 'TolX');
tol_abs = optimget(options, 'TolFun');
max_iter = optimget(options, 'MaxIter');

% initial condition
x = x0;
f0 = F(x0);
resnorm0 = norm(f0);
f = f0;
resnorm = resnorm0;

% Newton-Raphson looping
niter = 0;
while (resnorm>tol_abs && niter<max_iter)
    niter = niter + 1;
    
    J = dFdx(x);
    dx = - J \ f;
    x = x + dx;
    
    f = F(x);
    resnorm = norm(f);
end

retval = 0;
if (resnorm>tol_abs); retval = 1; end

return
end
