function DPCInitCap()
%DPCInitCap
% Solve cap center L0 = ASTART for inital pressure state
% Xinit = 3 * p_background
% so L0 must be the intersection of envelop and cap
% 

DPCGlobals;

cr = cr0;
xi = 3 * geop;

ok = 0;

if 0
    xo = 0;
    fo = cr * (aa-ac);
    x = xi + ab;
    for iter = 1:30
        f = x + cr*(aa-ac*exp(-ab*x)) - xi;
        if abs(f) < 1e-5
            ok = 1; break;
        end
        xn = x - f*(x-xo)/(f-fo);
        xo = x;
        fo = f;
        x = xn;
    end
else
    fun = @(ll) (ll-xi) + cr*(aa-ac*exp(-ab*ll));
    [x,~,ok] = fsolve(fun,0.0);
end

if ok == 1
    astart = x;
else
    error('astart failed to converge');
end

return
end

