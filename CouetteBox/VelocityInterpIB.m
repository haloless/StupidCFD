

function [u,v] = VelocityInterpIB(umac,vmac,nx,ny,dx,dy,dt)

EBGlobals;

u = umac;
v = vmac;

for i = 1:nx+3
for j = 1:ny+2
    if eb_uflag(i,j) == 1
        u(i,j) = eb_umac(i,j);
    elseif eb_uflag(i,j) == 2
        uint = 0.0;
        cnt = 0;
        if eb_uflag(i+1,j) == 0
            uint = uint + (umac(i+1,j)-eb_usurf(i+1,j)) / ebls_umac(i+1,j);
            cnt = cnt + 1;
        end
        if eb_uflag(i-1,j) == 0
            uint = uint + (umac(i-1,j)-eb_usurf(i-1,j)) / ebls_umac(i-1,j);
            cnt = cnt + 1;
        end
        if eb_uflag(i,j+1) == 0
            uint = uint + (umac(i,j+1)-eb_usurf(i,j+1)) / ebls_umac(i,j+1);
            cnt = cnt + 1;
        end
        if eb_uflag(i,j-1) == 0
            uint = uint + (umac(i,j-1)-eb_usurf(i,j-1)) / ebls_umac(i,j-1);
            cnt = cnt + 1;
        end
        if cnt > 0
            u(i,j) = eb_usurf(i,j) + ebls_umac(i,j) * (uint/cnt);
        else
            error(['Bad IB for U ',int2str(i),',',int2str(j)]);
        end
    end
end
end

for i = 1:nx+2
for j = 1:ny+3
    if eb_vflag(i,j) == 1
        v(i,j) = eb_vmac(i,j);
    elseif eb_vflag(i,j) == 2
        vint = 0.0;
        cnt = 0;
        if eb_vflag(i+1,j) == 0
            vint = vint + (vmac(i+1,j)-eb_vsurf(i+1,j)) / ebls_vmac(i+1,j);
            cnt = cnt + 1;
        end
        if eb_vflag(i-1,j) == 0
            vint = vint + (vmac(i-1,j)-eb_vsurf(i-1,j)) / ebls_vmac(i-1,j);
            cnt = cnt + 1;
        end
        if eb_vflag(i,j+1) == 0
            vint = vint + (vmac(i,j+1)-eb_vsurf(i,j+1)) / ebls_vmac(i,j+1);
            cnt = cnt + 1;
        end
        if eb_vflag(i,j-1) == 0
            vint = vint + (vmac(i,j-1)-eb_vsurf(i,j-1)) / ebls_vmac(i,j-1);
            cnt = cnt + 1;
        end
        if cnt > 0
            v(i,j) = eb_vsurf(i,j) + ebls_vmac(i,j) * (vint/cnt);
        else
            error(['Bad IB for V ',int2str(i),',',int2str(j)]);
        end
    end
end
end


return
end









