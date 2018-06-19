function [N] = fdpmShapeMLS(rx,ry,wgt)
%fdpmShapeMLS
% Calculate MLS shape function.
%
% Input:
% rx,ry: contains relative positions
% wgt: weights
%
% Output:
%


% support polynomial
P = fdpmShapeBasis(rx,ry);
% enhance by constant vector
P(:,2:end+1) = P;
P(:,1) = 1;

% center point
% if it is not the center point, then set [1, x,y, x^2/2,x*y,y^2/2]
ph = [ 1, 0,0, 0,0,0 ];

if 0
    % reduce to 1st order
    P = P(:,1:3);
    ph = [ 1, 0,0 ];
end

% moment matrix
Amat = P' * diag(wgt) * P;

%
acoef = (Amat' \ P') * diag(wgt);


% nodal shape function
N = ph * acoef;


return
end


function [P] = ShapeBasis(rx,ry)
%fdpmShapeBasis
% Return the polynomial basis
% Currently upto 2nd-order [x, y, x2, xy, y2]
%

rx = rx(:);
ry = ry(:);


% the order of polynomial
completeness = 2;

%
switch completeness
case {1}
	% x,y
	P = [ rx, ry ];
case {2}
	% x,y, x^2
	P = [ rx, ry, 1/2*rx.^2, rx.*ry, 1/2*ry.^2 ];
otherwise
	error('Unsupported completeness=%d',completeness);
end

return
end



function [w] = fdpmWeightFunc(q, cutoff)
%fdpmWeightFunc
% Calculate weight function 
%

if exist('cutoff','var')
	q = q ./ cutoff;
end

w = (q<1) .* (1 - 6*q.^2 + 8*q.^3 - 3*q.^4);

return
end






