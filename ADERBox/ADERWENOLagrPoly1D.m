
function [ ps ] = ADERWENOLagrPoly1D(qs,der)

ADERWENOGlobals1D;

% if (isrow(qs)); qs = qs'; end

d = 0;
if (nargin > 1); d = der; end

switch d
case 0
    ps = SumPoly(LagrPsi, qs);
case 1
    ps = SumPoly(LagrPsiDer1, qs);
otherwise
    error('Invalid derivative');
end

return
end


