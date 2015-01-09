

function [ Du Dv ] = EBVelocityDiffusion( ...
umac,vmac, ...
eb_umac,eb_vmac, ebls_umac, ebls_vmac, ...
nx,ny,dx,dy)


EBGlobals;


Du = zeros(nx+3,ny+2);
for i = 2:nx+2
for j = 2:ny+1
    uc = umac(i,j);
    if ebls_umac(i,j) > 0
        % XLO
        if ebls_umac(i-1,j) > 0
            hw = 1.0;
            uw = umac(i-1,j);
        else
            hw = ebls_umac(i,j) / (ebls_umac(i,j)-ebls_umac(i-1,j));
            uw = hw*eb_umac(i-1,j) + (1.0-hw)*eb_umac(i,j);
        end
        % XHI
        if ebls_umac(i+1,j) > 0
            he = 1.0;
            ue = umac(i+1,j);
        else
            he = ebls_umac(i,j) / (ebls_umac(i,j)-ebls_umac(i+1,j));
            ue = he*eb_umac(i+1,j) + (1.0-he)*eb_umac(i,j);
        end
        % YLO
        if ebls_umac(i,j-1) > 0
            hs = 1.0;
            us = umac(i,j-1);
        else
            hs = ebls_umac(i,j) / (ebls_umac(i,j)-ebls_umac(i,j-1));
            us = hs*eb_umac(i,j-1) + (1.0-hs)*eb_umac(i,j);
        end
        % YHI
        if ebls_umac(i,j+1) > 0
            hn = 1.0;
            un = umac(i,j+1);
        else
            hn = ebls_umac(i,j) / (ebls_umac(i,j)-ebls_umac(i,j+1));
            un = hn*eb_umac(i,j+1) + (1.0-hn)*eb_umac(i,j);
        end
        
        Du(i,j) = nu * ( ...
        1.0 / dx^2 * ((ue-uc)/he + (uw-uc)/hw) + ...
        1.0 / dy^2 * ((un-uc)/hn + (us-uc)/hs));
    end % (I,j) is fluid point
end
end

Dv = zeros(nx+2,ny+3);
for i = 2:nx+1
for j = 2:ny+2
    vc = vmac(i,j);
    if ebls_vmac(i,j) > 0
        % XLO
        if ebls_vmac(i-1,j) > 0
            hw = 1.0;
            vw = vmac(i-1,j);
        else
            hw = ebls_vmac(i,j) / (ebls_vmac(i,j)-ebls_vmac(i-1,j));
            vw = hw*eb_vmac(i-1,j) + (1.0-hw)*eb_vmac(i,j);
        end
        % XHI
        if ebls_vmac(i+1,j) > 0
            he = 1.0;
            ve = vmac(i+1,j);
        else
            he = ebls_vmac(i,j) / (ebls_vmac(i,j)-ebls_vmac(i+1,j));
            ve = he*eb_vmac(i+1,j) + (1.0-he)*eb_vmac(i,j);
        end
        % YLO
        if ebls_vmac(i,j-1) > 0
            hs = 1.0;
            vs = vmac(i,j-1);
        else
            hs = ebls_vmac(i,j) / (ebls_vmac(i,j)-ebls_vmac(i,j-1));
            vs = hs*eb_vmac(i,j-1) + (1.0-hs)*eb_vmac(i,j);
        end
        % YHI
        if ebls_vmac(i,j+1) > 0
            hn = 1.0;
            vn = vmac(i,j+1);
        else
            hn = ebls_vmac(i,j) / (ebls_vmac(i,j)-ebls_vmac(i,j+1));
            vn = hn*eb_vmac(i,j+1) + (1.0-hn)*eb_vmac(i,j);
        end
        
        Dv(i,j) = nu * ( ...
        1.0/dx^2 * ((ve-vc)/he + (vw-vc)/hw) + ...
        1.0/dy^2 * ((vn-vc)/hn + (vs-vc)/hs));
    end % (i,J) is fluid point
end
end

return
end


