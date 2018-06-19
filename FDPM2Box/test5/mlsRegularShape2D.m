function [N,Nx,Ny] = mlsRegularShape2D(xint, xs,ys, re, varargin)
%mlsRegularShape1D
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

%
if isrow(xs)
    xs = xs.';
end
if isrow(ys)
    ys = ys.';
end

rx = xint(1) - xs;
ry = xint(2) - ys;
% rr = sqrt(rx.^2 + ry.^2);


% set default values
mls_order = 1;
% mls_order = -1;
% mls_order = 2;
mls_volume = ones(size(xs));

mlsParseInterpOptions;


% weight and derivative
[w, wx,wy] = mlsRegularWeight2D(xint, xs,ys, re);
w = w .* mls_volume;
wx = wx .* mls_volume;
wy = wy .* mls_volume;

%
[P,Px,Py] = ShapeBasis(rx,ry, mls_order);
Pt = P.';
Pxt = Px.';
Pyt = Py.';

W = diag(w);

A = Pt * W * P;

invA = inv(A);

ah = invA * Pt * W;

[ph,phx] = ShapeBasis(0,0, mls_order);

N = ph * ah;

if nargout > 1
    % derivative
    
    Wx = diag(wx);
    Wy = diag(wy);
    
    Ax = Pxt * W * P + Pt * Wx * P + Pt * W * Px;
    Ay = Pyt * W * P + Pt * Wy * P + Pt * W * Py;
    
    invAx = -invA * Ax * invA;
    invAy = -invA * Ay * invA;
    
    ahx = invAx * Pt * W + invA * Pxt * W + invA * Pt * Wx;
    ahy = invAy * Pt * W + invA * Pyt * W + invA * Pt * Wy;
    
    Nx = ph * ahx;
    Ny = ph * ahy;
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [w,wx,wy] = WeightFunc(rx,ry,re);

% w = w .* mls_volume;
% wx = wx .* mls_volume;
% wy = wy .* mls_volume;

% [P,Px,Py] = ShapeBasis(rx,ry, mls_order);

% W = diag(w);
% Wx = diag(wx);
% Wy = diag(wy);

% A = P' * W * P;
% Ax = Px' * W * P + P' * Wx * P + P' * W * Px;
% Ay = Py' * W * P + P' * Wy * P + P' * W * Py;

% invA = inv(A);
% invAx = -invA * Ax * invA;
% invAy = -invA * Ay * invA;

% ah = invA * P' * W;
% ahx = invAx * P' * W + invA * Px' * W + invA * P' * Wx;
% ahy = invAy * P' * W + invA * Py' * W + invA * P' * Wy;

% % center point
% ph = ShapeBasis(0,0, mls_order);

% N = ph * ah;
% Nx = ph * ahx;
% Ny = ph * ahy;

% % TODO derivatives over 2nd-order are difficult to calculate!!!
% % They are here approximated
% if mls_order >= 2
    % Nxx = ah(4,:);
    % Nxy = ah(5,:);
    % Nyy = ah(6,:);
% else
    % Nxx = [];
    % Nxy = [];
    % Nyy = [];
% end



return
end


function [P,Px,Py] = ShapeBasis(rx, ry, completeness)
%ShapeBasis
% Return the polynomial basis
% Currently upto 2nd-order [x, y, x2, xy, y2]
%

vec0 = zeros(size(rx));
vec1 = ones(size(rx));


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
    P = [ vec1, rx, ry ];
    Px = [ vec0, vec1, vec0 ];
    Py = [ vec0, vec0, vec1 ];
case {-1}
    % 1,x,y,xy
    P = [ vec1, rx, ry, rx.*ry ];
    Px = [ vec0, vec1, vec0, ry ];
    Py = [ vec0, vec0, vec1, rx ];
case {2}
    % 1, x,y, x^2, x*y, y^2
    P = [ vec1, rx, ry, 1/2*rx.^2, rx.*ry, 1/2*ry.^2 ];
    Px = [ vec0, vec1, vec0, rx, ry, vec0 ];
    Py = [ vec0, vec0, vec1, vec0, rx, ry ];
otherwise
	error('Unsupported completeness=%d',completeness);
end

return
end







