function [mob] = EwaldMobMatrix(ewald,np,xp)

xi = ewald.xi;
xi2 = xi^2;
lx = ewald.lx;
ly = ewald.ly;
lz = ewald.lz;

% currently FT version
ver = 6;

mob = zeros(np*ver,np*ver);


%
% self part
%
for i = 1:np
    mob = AddMobMatrix(mob, np, i,i, 0,0,0, ...
    ewald.self_a,ewald.self_a, 0.0, ewald.self_c,ewald.self_c, ...
    0.0,0.0, 0.0, 0.0,0.0,0.0);
end

% save 
mdiag = mob;

%
% Ewald Sum in real space
%
for i = 1:np
for j = 1:np
    
    xdiff = xp(:,j) - xp(:,i);
    
    for m1 = -ewald.rnumx:ewald.rnumx
    for m2 = -ewald.rnumy:ewald.rnumy
    for m3 = -ewald.rnumz:ewald.rnumz
        rlx = ewald.lx * m1;
        rly = ewald.ly * m2;
        rlz = ewald.lz * m3;
        
        xx = xdiff(1) + rlx;
        yy = xdiff(2) + rly;
        zz = xdiff(3) + rlz;
        rr = sqrt(xx^2 + yy^2 + zz^2);
        
        if (rr == 0) % exclude self part
            continue;
        elseif (rr > ewald.rmax)
            continue;
        end
        
        ex = xx / rr;
        ey = yy / rr;
        ez = zz / rr;
        
        [xa,ya,yb,xc,yc] = EwaldScalarReal(xi,rr);
        
        mob = AddMobMatrix(mob,np,i,j,ex,ey,ez, ...
        xa,ya, yb, xc,yc, 0,0,0,0,0,0);
    end
    end
    end
end
end

% save
mreal = mob;


%
% Ewald Sum in wave space
%
k0x = 2*pi / ewald.lx;
k0y = 2*pi / ewald.ly;
k0z = 2*pi / ewald.lz;
for m1 = -ewald.knumx:ewald.knumx
for m2 = -ewald.knumy:ewald.knumy
for m3 = -ewald.knumz:ewald.knumz
    k1 = k0x * m1;
    k2 = k0x * m2;
    k3 = k0x * m3;
    kk = k1*k1 + k2*k2 + k3*k3;
    k = sqrt(kk);
    
    if (kk == 0)
        continue;
    elseif (k > ewald.kmax)
        continue;
    end
    
    k4z = kk / 4 / xi2;
    kexp = ewald.pivol * (1.0 + k4z * (1.0 + 2.0 * k4z)) / kk * exp(- k4z);
    
    ex = k1 / k;
    ey = k2 / k;
    ez = k3 / k;
    
    for i = 1:np
    for j = 1:np
        ya = 6.0 * (1.0 - kk / 3.0) * kexp;
        yb = 3.0 * k * kexp;
        yc = 3.0 / 2.0 * kk * kexp;
        yg = 3.0 * (1.0 - 4.0 / 15.0 * kk) * k * kexp;
        yh = 3.0 / 2.0 * kk * kexp;
        ym = 3.0 * (1.0 - kk / 5.0) * kk * kexp;
        
        xx = xp(1,j) - xp(1,i);
        yy = xp(2,j) - xp(2,i);
        zz = xp(3,j) - xp(3,i);
        
        cf = cos(k1*xx + k2*yy + k3*zz);
        sf = -sin(k1*xx + k2*yy + k3*zz);
        
        mob = AddMobMatrix(mob,np,i,j,ex,ey,ez, ...
        0.0,cf*ya, sf*yb, 0.0,cf*yc, ...
        0,0,0,0,0,0);
    end
    end
end
end
end



return
end


