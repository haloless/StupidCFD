
clear all;

I2 = eye(2);

% xc = 2.0;
% yc = 0.0;
xc = 2.33;
yc = 1.37;
rad = 1.0;

func_phi = @(x,y) sqrt(x^2+y^2);
func_phix = @(x,y) - x / sqrt(x^2+y^2);
func_phiy = @(x,y) - y / sqrt(x^2+y^2);

if (1)
    % calculate integral
    nseg = 240;
    dl = 2*pi*rad / nseg;
    
    fint = zeros(2,1);
    for i = 1:nseg
        theta = 2*pi/nseg * (i-1);
        x = xc + rad*cos(theta);
        y = yc + rad*sin(theta);
        
        Evec = [ func_phix(x,y); func_phiy(x,y) ];
        
        sigma = Evec*Evec' - 0.5*dot(Evec,Evec)*I2;
        
        nvec = [ cos(theta); sin(theta) ];
        
        fvec = sigma * nvec;
        
        fint = fint + fvec.*dl;
    end
    
    fint
end

xlo = xc - rad*2;
xhi = xc + rad*2;
ylo = yc - rad*2;
yhi = yc + rad*2;

nn = 16;
% nn = 32;
% nn = 64;
% nn = 128;
% nn = 256;
% nn = 512;
nx = nn;
ny = nn;
dx = (xhi-xlo) / nx;
dy = (yhi-ylo) / ny;
dv = dx * dy;
dsx = dv / dx;
dsy = dv / dy;

xcell = zeros(nx,ny);
ycell = zeros(nx,ny);
phi = zeros(nx,ny);
sdf = zeros(nx,ny);
for j = 1:ny
for i = 1:nx
    x = xlo + (i-0.5)*dx;
    y = ylo + (j-0.5)*dy;
    rc = sqrt((x-xc)^2+(y-yc)^2);
    
    xcell(i,j) = x;
    ycell(i,j) = y;
    sdf(i,j) = rc - rad;
    if sdf(i,j) > 0
        phi(i,j) = func_phi(x,y);
    else
        phi(i,j) = 0.0;
    end
end
end

if (1)
    figure;
    contourf(xcell,ycell,phi);
    axis equal;
    axis([xlo, xhi, ylo, yhi]);
    % figure;
    % contourf(xcell,ycell,sdf);
end

fsumx = 0;
fsumy = 0;

if (1)
    % calculate -grad(phi), use 2nd-order FD
    gphix = zeros(nx,ny);
    gphiy = zeros(nx,ny);
    
    for j = 2:ny-1
    for i = 2:nx-1
        if sdf(i,j) > 0.0
            if sdf(i-1,j) <= 0
                phix = -(-3*phi(i,j)+4*phi(i+1,j)-phi(i+2,j)) / (dx*2);
            elseif sdf(i+1,j) <= 0
                phix = (-3*phi(i,j)+4*phi(i-1,j)-phi(i-2,j)) / (dx*2);
            else
                phix = -(phi(i+1,j)-phi(i-1,j)) / (dx*2);
            end
            
            if sdf(i,j-1) <= 0
                phiy = -(-3*phi(i,j)+4*phi(i,j+1)-phi(i,j+2)) / (dy*2);
            elseif sdf(i,j+1) <= 0
                phiy = (-3*phi(i,j)+4*phi(i,j-1)-phi(i,j-2)) / (dy*2);
            else
                phiy = -(phi(i,j+1)-phi(i,j-1)) / (dy*2);
            end
            
            gphix(i,j) = phix;
            gphiy(i,j) = phiy;
        end
    end
    end
end

if (0)
    for j = 2:ny-1
    for i = 2:nx-1
        if sdf(i,j) > 0.0
            okx = 1;
            oky = 1;
            
            if sdf(i-1,j) <= 0
                % phix = -(phi(i+1,j)-phi(i,j)) / dx;
                phix = -(-3*phi(i,j)+4*phi(i+1,j)-phi(i+2,j)) / (dx*2);
                nvecx = 1;
            elseif sdf(i+1,j) <= 0
                % phix = -(phi(i,j)-phi(i-1,j)) / dx;
                phix = (-3*phi(i,j)+4*phi(i-1,j)-phi(i-2,j)) / (dx*2);
                nvecx = -1;
            else
                phix = -(phi(i+1,j)-phi(i-1,j)) / (dx*2);
                nvecx = 0;
                okx = 0;
            end
            
            if sdf(i,j-1) <= 0
                % phiy = -(phi(i,j+1)-phi(i,j)) / dy;
                phiy = -(-3*phi(i,j)+4*phi(i,j+1)-phi(i,j+2)) / (dy*2);
                nvecy = 1;
            elseif sdf(i,j+1) <= 0
                % phiy = -(phi(i,j)-phi(i,j-1)) / dy;
                phiy = (-3*phi(i,j)+4*phi(i,j-1)-phi(i,j-2)) / (dy*2);
                nvecy = -1;
            else
                phiy = -(phi(i,j+1)-phi(i,j-1)) / (dy*2);
                nvecy = 0;
                oky = 0;
            end
            
            if (1)
                % set to analytical
                x = xcell(i,j);
                y = ycell(i,j);
                phix = func_phix(x,y);
                phiy = func_phiy(x,y);
            end
            Evec = [ phix; phiy ];
            sigma = Evec*Evec' - 0.5*dot(Evec,Evec)*I2;
            
            if okx
                fsumx = fsumx + sigma(1,1)*nvecx*dsx;
                fsumy = fsumy + sigma(2,1)*nvecx*dsx;
            end
            
            if oky
                fsumx = fsumx + sigma(1,2)*nvecy*dsy;
                fsumy = fsumy + sigma(2,2)*nvecy*dsy;
            end
        end
    end
    end
end

if (0)
    % use height fraction, analytical -grad(phi)
    for j = 2:ny-1
    for i = 2:nx-1
        if sdf(i,j) > 0
            okx = 0;
            if sdf(i-1,j) <= 0
                okx = 1;
                hfrac = sdf(i,j) / (sdf(i,j)-sdf(i-1,j));
                xx = xcell(i,j) - hfrac*dx;
                yy = ycell(i,j);
                nvecx = 1;
            elseif sdf(i+1,j) <= 0
                okx = 1;
                hfrac = sdf(i,j) / (sdf(i,j)-sdf(i+1,j));
                xx = xcell(i,j) + hfrac*dx;
                yy = ycell(i,j);
                nvecx = -1;
            end
            if okx
                phix = func_phix(xx,yy);
                phiy = func_phiy(xx,yy);
                Evec = [ phix; phiy ];
                sigma = Evec*Evec' - 0.5*dot(Evec,Evec)*I2;
                fsumx = fsumx + sigma(1,1)*nvecx*dsx;
                fsumy = fsumy + sigma(2,1)*nvecx*dsx;
            end
            
            oky = 0;
            if sdf(i,j-1) <= 0
                oky = 1;
                hfrac = sdf(i,j) / (sdf(i,j)-sdf(i,j-1));
                xx = xcell(i,j);
                yy = ycell(i,j) - hfrac*dy;
                nvecy = 1;
            elseif sdf(i,j+1) <= 0
                oky = 1;
                hfrac = sdf(i,j) / (sdf(i,j)-sdf(i,j+1));
                xx = xcell(i,j);
                yy = ycell(i,j) + hfrac*dy;
                nvecy = -1;
            end
            if oky
                phix = func_phix(xx,yy);
                phiy = func_phiy(xx,yy);
                Evec = [ phix; phiy ];
                sigma = Evec*Evec' - 0.5*dot(Evec,Evec)*I2;
                fsumx = fsumx + sigma(1,2)*nvecy*dsy;
                fsumy = fsumy + sigma(2,2)*nvecy*dsy;
            end
        end
    end
    end
end

if (1)
    % use height fraction, directional extrapolation
    for j = 2:ny-1
    for i = 2:nx-1
        if sdf(i,j) > 0
            okx = 0;
            if sdf(i-1,j) <= 0
                okx = 1;
                hfrac = sdf(i,j) / (sdf(i,j)-sdf(i-1,j));
                phix = (1+hfrac)*gphix(i,j) - hfrac*gphix(i+1,j);
                phiy = (1+hfrac)*gphiy(i,j) - hfrac*gphiy(i+1,j);
                nvecx = 1;
            elseif sdf(i+1,j) <= 0
                okx = 1;
                hfrac = sdf(i,j) / (sdf(i,j)-sdf(i+1,j));
                phix = (1+hfrac)*gphix(i,j) - hfrac*gphix(i-1,j);
                phiy = (1+hfrac)*gphiy(i,j) - hfrac*gphiy(i-1,j);
                nvecx = -1;
            end
            if okx
                Evec = [ phix; phiy ];
                sigma = Evec*Evec' - 0.5*dot(Evec,Evec)*I2;
                fsumx = fsumx + sigma(1,1)*nvecx*dsx;
                fsumy = fsumy + sigma(2,1)*nvecx*dsx;
            end
            
            oky = 0;
            if sdf(i,j-1) <= 0
                oky = 1;
                hfrac = sdf(i,j) / (sdf(i,j)-sdf(i,j-1));
                phix = (1+hfrac)*gphix(i,j) - hfrac*gphix(i,j+1);
                phiy = (1+hfrac)*gphiy(i,j) - hfrac*gphiy(i,j+1);
                nvecy = 1;
            elseif sdf(i,j+1) <= 0
                oky = 1;
                hfrac = sdf(i,j) / (sdf(i,j)-sdf(i,j+1));
                phix = (1+hfrac)*gphix(i,j) - hfrac*gphix(i,j-1);
                phiy = (1+hfrac)*gphiy(i,j) - hfrac*gphiy(i,j-1);
                nvecy = -1;
            end
            if oky
                Evec = [ phix; phiy ];
                sigma = Evec*Evec' - 0.5*dot(Evec,Evec)*I2;
                fsumx = fsumx + sigma(1,2)*nvecy*dsy;
                fsumy = fsumy + sigma(2,2)*nvecy*dsy;
            end
        end
    end
    end
end



disp(['R/dx=',num2str(rad/dx),'; D/dx=',num2str(rad*2/dx)]);
[fsumx,fsumy]




