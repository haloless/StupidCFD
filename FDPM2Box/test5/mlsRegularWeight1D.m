function [w, wx] = mlsRegularWeight1D(xint, xs, re)

%
ee = 1.0e-3;
% ee = 1.0e-4;
% ee = 1.0e-5;
gg = 2;

assert(gg>=2 && mod(gg,2)==0);

e1 = (1+ee)^(-2);
e2 = ee^(-2) - e1;



nint = size(xint, 1);
npnt = length(xs(:));

for i = 1:nint
    xx = xint(i);
    
    dd = xx - xs;
    rr = abs(dd);
    ss = rr ./ re;
    
    % neighbors
    ok = ss < 1;
    mm = find(ok);
    nn = length(mm);
    
    %
    sg = zeros(size(xs));
    sg(mm) = ss(mm).^gg + ee;
    
    what = zeros(size(xs));
    what(mm) = (sg(mm).^(-2) - e1) ./ e2;
    
    wall = sum(what(mm));
    w = what ./ wall;
    
    
    %
    % compute derivative
    %
    if nargout > 1
        
        % s * ds/dx = (x-xk)/re^2
        sdsdx = dd(mm) ./ re^2;
        
        whatx = zeros(size(xs));
        whatx(mm) = -2.0/e2 .* sg(mm).^(-3) * gg .* ss(mm).^(gg-2) .* sdsdx;
        
        wallx = sum(whatx(mm));
        
        wx = zeros(size(xs));
        wx(mm) = (whatx(mm).*wall - what(mm).*wallx) ./ wall^2;
    end
    
    
end

return
end

