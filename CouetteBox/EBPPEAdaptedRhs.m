

function [ rhs ] = EBPPEAdaptedRhs( ...
ustar,vstar, ...
ebls, ebls_umac, ebls_vmac, eb_umac, eb_vmac, ...
nx,ny,dx,dy,dt)

EBGlobals;

divu = zeros(nx+2,ny+2);

for i = 2:nx+1
for j = 2:ny+1
    if ebls(i,j) > 0
        % XLO
        if ebls(i-1,j) > 0
            xl = 0;
            ul = ustar(i,j);
        else
            if ebls_umac(i,j) > 0 % extend
                if ebls_umac(i-1,j)>0; error('EB invalid'); end
                theta = ebls_umac(i,j) / (ebls_umac(i,j)-ebls_umac(i-1,j));
                xl = -theta;
                ul = theta*eb_umac(i-1,j) + (1.0-theta)*eb_umac(i,j);
            else % shrink
                if ebls_umac(i+1,j)<=0; error('EB invalid'); end
                theta = ebls_umac(i+1,j) / (ebls_umac(i+1,j)-ebls_umac(i,j));
                xl = 1.0 - theta;
                ul = theta*eb_umac(i,j) + (1.0-theta)*eb_umac(i+1,j);
            end
        end
        % XHI
        if ebls(i+1,j) > 0
            xr = 0;
            ur = ustar(i+1,j);
        else
            if ebls_umac(i+1,j) > 0 % extend
                if ebls_umac(i+2,j)>0; error('EB invalid'); end
                theta = ebls_umac(i+1,j) / (ebls_umac(i+1,j)-ebls_umac(i+2,j));
                xr = theta;
                ur = theta*eb_umac(i+2,j) + (1.0-theta)*eb_umac(i+1,j);
            else % shrink
                if ebls_umac(i,j)<=0; error('EB invalid'); end
                theta = ebls_umac(i,j) / (ebls_umac(i,j)-ebls_umac(i+1,j));
                xr = -(1.0 - theta);
                ur = theta*eb_umac(i+1,j) + (1.0-theta)*eb_umac(i,j);
            end
        end
        % YLO
        if ebls(i,j-1) > 0
            yl = 0;
            vl = vstar(i,j);
        else
            if ebls_vmac(i,j) > 0 % extend
                if ebls_vmac(i,j-1) > 0; error('EB invalid'); end
                theta = ebls_vmac(i,j) / (ebls_vmac(i,j)-ebls_vmac(i,j-1));
                yl = -theta;
                vl = theta*eb_vmac(i,j-1) + (1.0-theta)*eb_vmac(i,j);
            else % shrink
                if ebls_vmac(i,j+1)<=0; error('EB invalid'); end
                theta = ebls_vmac(i,j+1) / (ebls_vmac(i,j+1)-ebls_vmac(i,j));
                yl = 1.0 - theta;
                vl = theta*eb_vmac(i,j) + (1.0-theta)*eb_vmac(i,j+1);
            end
        end
        % YHI
        if ebls(i,j+1) > 0
            yr = 0;
            vr = vstar(i,j+1);
        else
            if ebls_vmac(i,j+1) > 0 % extend
                if ebls_vmac(i,j+2)>0; error('EB invalid'); end
                theta = ebls_vmac(i,j+1) / (ebls_vmac(i,j+1)-ebls_vmac(i,j+2));
                yr = theta;
                vr = theta*eb_vmac(i,j+2) + (1.0-theta)*eb_vmac(i,j+1);
            else % shrink
                if ebls_vmac(i,j)<=0; error('EB invalid'); end
                theta = ebls_vmac(i,j) / (ebls_vmac(i,j)-ebls_vmac(i,j+1));
                yr = -(1.0 - theta);
                vr = theta*eb_vmac(i,j+1) + (1.0-theta)*eb_vmac(i,j);
            end
        end
        
        % xl = 0; ul = ustar(i,j);
        % xr = 0; ur = ustar(i+1,j);
        % yl = 0; vl = vstar(i,j);
        % yr = 0; vr = vstar(i,j+1);
        
        hx = xr - xl + 1.0;
        hy = yr - yl + 1.0;
        % hx = 1.0;
        % hy = 1.0;
        % beta = min([1.0, hx, hy]);
        beta = 1;
        divu(i,j) = ((ur-ul)/(hx*dx) + (vr-vl)/(hy*dy)) * beta;
    end % (i,j) is fluid cell
end
end

rhs = divu(2:nx+1,2:ny+1);
rhs = -rho/dt .* rhs;
rhs = reshape(rhs, nx*ny,1);
rhs = rhs - sum(rhs)/(nx*ny);

return
end



