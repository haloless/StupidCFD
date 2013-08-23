% ## Copyright (C) 2013 admin_2
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

% ## MLS_shape

% ## Author: admin_2 <admin_2@KOSHIZUKA>
% ## Created: 2013-07-19

function [ phi,phi_x,phi_y ] = MLS_shape (xc,yc,xs,ys,re,neigh)
% Description:
% return phi, dphi_dx, dphi_dy in sparse row vectors

use_pseudo_inverse = 1;

% weight function, in sparse vector
[w,dwdx,dwdy] = MLS_weight(xc,yc,xs,ys,re,neigh);
% convert to sparse matrix
W = diag(w);
W_x = diag(dwdx);
W_y = diag(dwdy);

% basis for (xc,yc)
[pc,pc_x,pc_y] = MLS_basis(xc,yc);
pc = sparse(pc);
pc_x = sparse(pc_x);
pc_y = sparse(pc_y);

% all known basis for (xi,yi)
P = MLS_basis(xs,ys,neigh)';

% building blocks for MLS shape func.
B = P' * W;
A = B * P;
if (use_pseudo_inverse)
    % use pseudo-inverse instead
    [U,S,V] = svd(A);
    smax = max(abs(diag(S)));
    trunc = find(abs(diag(S))<smax*1e-6);
    invS = 1 ./ diag(S);
    invS(trunc) = 0;
    invS = diag(invS);
    invA = V * invS * U';
else
    % no-good if ill-conditioned
    invA = inv(A);
end

B_x = P' * W_x;
B_y = P' * W_y;
A_x = B_x * P;
A_y = B_y * P;
invA_x = -invA * A_x * invA;
invA_y = -invA * A_y * invA;

% MLS shape func., in sparse row vector
phi = pc' * invA * B;

% MLS shape func. derivatives
phi_x = pc_x'*invA*B + pc'*invA_x*B + pc'*invA*B_x;
phi_y = pc_y'*invA*B + pc'*invA_y*B + pc'*invA*B_y;


return
end



