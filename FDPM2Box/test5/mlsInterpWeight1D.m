function [w, wx] = mlsInterpWeight1D(xint, xs, re)

% the power index
% a = 2;
a = 4;


nint = size(xint, 1);
npnt = length(xs(:));

for i = 1:nint
    xx = xint(i);
    
    dd = xx - xs;
    rr = abs(dd);
    
    ss = rr ./ re;
    [smin, imin] = min(ss);
    smina = smin.^a;
    
    % neighbors
    ok = ss < 1;
    mm = length(find(ok));
    
    sna = ss.^(-a);
    sna(imin) = 0;
    sna(~ok) = 0;
    
    sa = smina .* sna;
    sa(imin) = 1;
    sa(~ok) = 0;
    
    sall = sum(sa) + smina*mm;
    
    w = (sa - smina) ./ sall;
    w(~ok) = 0;
    
    %
    % compute derivative
    %
    if nargout > 1
        
        dsdx = 1/re .* (dd ./ rr);
        if smin == 0
            dsdx(imin) = 0;
        end
        dsdx(~ok) = 0;
        
        dsall = a * smin^(a-1) * dsdx(imin) * (sum(sna)+mm);
        for j = 1:npnt
            if ok(j) && j~=imin
                dsall = dsall - smin^a * a * ss(j)^(-a-1) * dsdx(j);
            end
        end
        
        wx = zeros(size(w));
        for j = 1:npnt
            if ok(j)
                if j == imin
                    wx(j) = -a * ss(j)^(a-1) * dsdx(j)*sall - (1-ss(j)^a)*dsall;
                    wx(j) = wx(j) / sall^2;
                else
                    wx(j) = a*ss(imin)^(a-1)*dsdx(imin)*(ss(j)^(-a)-1);
                    wx(j) = wx(j) + ss(imin)^a * (-a) * ss(j)^(-a-1) * dsdx(j);
                    wx(j) = wx(j)*sall - ss(imin)^a*(ss(j)^(-a)-1)*dsall;
                    wx(j) = wx(j) / sall^2;
                end
            end
        end
        % sall
        % dsall
    end
    
    
end

return
end

