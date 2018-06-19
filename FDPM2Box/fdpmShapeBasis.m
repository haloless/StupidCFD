function [P] = fdpmShapeBasis(rx,ry)
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

