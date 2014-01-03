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
% ## Created: 2013-09-04

function [ ret,sol ] = EBPPESolver_cg (nx,ny,dx,dy, adiag,bx,by, ...
sol,rhs, eps_rel,eps_abs,max_iter)
% Description
% EBGlobals;
SMALL_NUMBER = eps * 16;
%
STAT_GOOD = 0;
STAT_LOSS_ACCURACY = 1;
STAT_ITER_UNSTAB = 2;
STAT_ITER_EXCEED = 8;
%
unstab_crit = 10;

% generate PPE coefficients
[be,bw,bn,bs,bp] = EBPPEGenerateCoef(nx,ny,dx,dy, adiag,bx,by);

%
% precond = 0;
precond = 1;
switch (precond)
case {1}
    % M = bp;
    % M(find(abs(M)<=eps)) = 1;
    % P = zeros(nx+2,ny+2);
    % P(2:nx+1,2:ny+1) = 1 ./ M;
    M = zeros(nx+2,ny+2);
    M(2:nx+1,2:ny+1) = 1 ./ bp;
otherwise
% no preconditioning
end


% CG solving
ret = STAT_GOOD;
conv = 0;
verbose = 1;
% sorig = zeros(nx+2,ny+2);
% sorig(2:nx+1,2:ny+1) = sol;
% sol(:,:) = 0;
x = zeros(nx+2,ny+2);
x(2:nx+1,2:ny+1) = sol;
%
b = zeros(nx+2,ny+2);
b(2:nx+1,2:ny+1) = rhs;
%
% r = b - applyOp(nx,ny, adiag,be,bw,bn,bs,bp, x);
r = b - EBPPEApplyOp(nx,ny, be,bw,bn,bs,bp, x);
rnorm = norm_inf(nx,ny, r);
rnorm0 = rnorm;
rnorm_min = rnorm;
rh_norm = rnorm0;

% initial check
if (rnorm==0.0 || rnorm<eps_rel*rh_norm || rnorm<eps_abs)
    ret = STAT_GOOD;
    return;
end

% set to residual-correction-form
x(2:nx+1,2:ny+1) = 0;

for nit = 1:max_iter
    % do precondition here
    switch (precond)
    case {1} % Jacobi
        % z = P .* r;
        z = M .* r;
    otherwise % no precondition
        z = r;
    end
    rho = dotxy(nx,ny, r,z);
    
    if (nit == 1)
        p = z;
    else
        beta = rho / rho_old;
        p = z + beta*p;
    end
    
    % q = applyOp(nx,ny, adiag,be,bw,bn,bs,bp, p);
    q = EBPPEApplyOp(nx,ny, be,bw,bn,bs,bp, p);
    pw = dotxy(nx,ny, p,q);
    if (pw==0.0)
        % we have a zero product, which means loss-of-accuracy
        ret = STAT_LOSS_ACCURACY;
        break;
    end
    
    alpha = rho / pw;
    x = x + alpha*p;
    % if (mod(nit,25)==0)
        % % r = b - applyOp(nx,ny, adiag,be,bw,bn,bs,bp, x);
        % r = b - EBPPEApplyOp(nx,ny, be,bw,bn,bs,bp, x);
    % else
        % r = r - alpha*q;
    % end
    r = r - alpha*q;
    rnorm = norm_inf(nx,ny,r);
    
    % convergence
    if (rnorm<eps_rel*rnorm0 || rnorm<eps_abs)
        conv = 1;
        break;
    end
    % stability
    if (rnorm>unstab_crit*rnorm_min)
        ret = STAT_ITER_UNSTAB;
        break;
    elseif (rnorm < rnorm_min)
        rnorm_min = rnorm;
    end
    
    %
    rho_old = rho;
end % end of CG loop

if (ret==STAT_GOOD && ~conv)
    warning('CG solver: failed to converge!');
    ret = STAT_ITER_EXCEED;
end

if (ret==STAT_GOOD || ret==STAT_ITER_EXCEED)
    sol = sol + x(2:nx+1,2:ny+1);
end
% sol = x(2:nx+1,2:ny+1);
return
end
% end of EBPPESolver

% function [ ret ] = applyBC(nx,ny, phi)
% ret = PressureBC(phi,nx,ny);
% return
% end

% function [ ret ] = applyOp(nx,ny, A,Be,Bw,Bn,Bs,Bp, phi)
% phi = applyBC(nx,ny, phi);
% I = 2:nx+1;
% J = 2:ny+1;
% ret = zeros(nx+2,ny+2);
% % ret(I,J) = A.*phi(I,J) + Bp.*phi(I,J) ...
% % - Be.*phi(I+1,J) - Bw.*phi(I-1,J) - Bn.*phi(I,J+1) - Bs.*phi(I,J-1);
% ret(I,J) = A.*phi(I,J) ...
% + Be.*(phi(I,J)-phi(I+1,J)) + Bw.*(phi(I,J)-phi(I-1,J)) ...
% + Bn.*(phi(I,J)-phi(I,J+1)) + Bs.*(phi(I,J)-phi(I,J-1));
% return
% end

function [ ret ] = dotxy(nx,ny, x,y)
I = 2:nx+1;
J = 2:ny+1;
ret = sum(sum(x(I,J).*y(I,J)));
return
end

function [ ret ] = norm_inf(nx,ny, x)
ret = norm(reshape(x(2:nx+1,2:ny+1),[],1),Inf);
return
end
