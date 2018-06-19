function [N,Nx,Ny, Nxx,Nxy,Nyy] = fdpmShapeMLS3(xc,yc,xs,ys,re)
%fdpmShapeMLS3
% Calculate MLS shape function.
%
% Input:
% rx,ry: contains relative positions (xj - x, yj - y)
% 
% wgt: weights
%
% Output:
%

% note the minus sign
% (xj-x) -> (x-xj)
% rx = -rx(:);
% ry = -ry(:);

mls_order = 1;
% mls_order = 2;

[w,wx,wy] = WeightFunc(xc-xs,yc-ys,re);

[P,~,~] = ShapeBasis(xs,ys, mls_order);

W = diag(w);
Wx = diag(wx);
Wy = diag(wy);

A = P' * W * P;
% Ax = Px' * W * P + P' * Wx * P + P' * W * Px;
% Ay = Py' * W * P + P' * Wy * P + P' * W * Py;
Ax = P' * Wx * P;
Ay = P' * Wy * P;

invA = inv(A);
invAx = -invA * Ax * invA;
invAy = -invA * Ay * invA;

B = P' * W;
Bx = P' * Wx;
By = P' * Wy;

ah = invA * B;
ahx = invAx * B + invA * Bx;
ahy = invAy * B + invA * By;

% center point
[ph,phx,phy] = ShapeBasis(xc,yc, mls_order);

N = ph * ah;
Nx = ph * ahx + phx*ah;
Ny = ph * ahy + phy*ah;

% TODO derivatives over 2nd-order are difficult to calculate!!!
% They are here approximated
if mls_order >= 2
    Nxx = ah(4,:);
    Nxy = ah(5,:);
    Nyy = ah(6,:);
else
    Nxx = [];
    Nxy = [];
    Nyy = [];
end



return
end


function [P,Px,Py] = ShapeBasis(x,y, completeness)
%ShapeBasis
% Return the polynomial basis
% Currently upto 2nd-order [x, y, x2, xy, y2]
%

vec0 = zeros(size(x));
vec1 = ones(size(x));


% the order of polynomial
% completeness = 1;
% completeness = 2;

%
switch completeness
case {0}
    % 1
    P = [ vec1 ];
    Px = [ vec0 ];
    Py = [ vec0 ];
case {1}
    % 1,x,y
    P = [ vec1, x, y ];
    Px = [ vec0, vec1, vec0 ];
    Py = [ vec0, vec0, vec1 ];
% case {2}
    % % x,y, x^2
    % P = [ vec1, x, y, 1/2*x.^2, x.*y, 1/2*y.^2 ];
    % Px = [ vec0, vec1, vec0, x, y, vec0 ];
    % Py = [ vec0, vec0, vec1, vec0, x, y ];
otherwise
	error('Unsupported completeness=%d',completeness);
end

% Px = zeros(size(P));
% Py = zeros(size(P));

return
end



function [w,wx,wy] = WeightFunc(rx,ry, cutoff)
%fdpmWeightFunc
% Calculate weight function 
%

rsmall = cutoff * 1.0e-6;

r = sqrt(rx.^2 + ry.^2);

q = r ./ cutoff;

ok = (q<1);

% w = w(q)
w = ok .* (1 - 6*q.^2 + 8*q.^3 - 3*q.^4);

if nargout > 1
    % dw/dq
    wq = ok .* (-12*q + 24*q.^2 - 12*q.^3);
    
    % dw/dr
    wr = wq ./ cutoff;
    
    % dw/dx, dw/dy
    wx = rx ./ (r+rsmall) .* wr;
    wy = ry ./ (r+rsmall) .* wr;
end

return
end






