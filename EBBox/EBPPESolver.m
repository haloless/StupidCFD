% ## Copyright (C) 2013 homu
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

% ## EBPPESolver

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-09-06

function [ ret,sol ] = EBPPESolver (nx,ny,dx,dy, adiag,bx,by, ...
sol,rhs, method,params)

% generate PPE coefficients
[Be,Bw,Bn,Bs,Bp] = EBPPEGenerateCoef(nx,ny,dx,dy, adiag,bx,by);

%
x = reshape(sol, nx*ny,1);
b = reshape(rhs, nx*ny,1);
%
xbuf = zeros(nx+2,ny+2);
afunc = @A_op;
%
switch (params.precond)
case {'jacobi'}
    invD = Bp;
    invD(abs(invD)<eps) = 1;
    invD = reshape(1 ./ invD, nx*ny,1);
    mfunc = @M_Jacobi_op;
otherwise % no preconditioner
    mfunc = [];
end

%
tol = params.tol;
maxit = params.maxit;
switch (method)
case {'pcg'}
    [x,ret,res,iter] = pcg(afunc,b, tol,maxit, mfunc);
case {'bicgstab'}
    [x,ret,res,iter] = bicgstab(afunc,b, tol,maxit, mfunc);
case {'gmres'} % OMG this GMRES is terribly slow...
    restart = params.restart;
    [x,ret,res,iter] = gmres(afunc,b, restart,tol,maxit, mfunc);
otherwise
    error('Unknown solving method: %s', method);
end

% result
sol = reshape(x, nx,ny);


% function handlers as Matrix operators
function [ y ] = A_op(x)
    xbuf(2:nx+1,2:ny+1) = reshape(x,nx,ny);
    xbuf = EBPPEApplyOp(nx,ny, Be,Bw,Bn,Bs,Bp, xbuf);
    y = reshape(xbuf(2:nx+1,2:ny+1), nx*ny,1);
end

% preconditioner
%
function [ z ] = M_Jacobi_op(r)
    z = invD .* r;
end
return
end


