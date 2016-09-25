
clear all;


R = 1.0;

xlo = -2.0;
xhi = 2.0;
ylo = -2.0;
yhi = 2.0;

nx = 40;
ny = 40;
dx = (xhi-xlo) / nx;
dy = (yhi-ylo) / ny;

xnode = linspace(xlo,xhi,nx+1);
ynode = linspace(ylo,yhi,ny+1);
[xnode,ynode] = ndgrid(xnode,ynode);

unode = zeros(nx+1,ny+1);
for j = 1:ny+1
for i = 1:nx+1
    xx = xnode(i,j);
    yy = ynode(i,j);
    rr = sqrt(xx^2+yy^2);
    
    if abs(rr-R) < 0.5*R
        nvecx = xx / rr;
        nvecy = yy / rr;
        dist = R - rr;
        % unode(i,j) = CubeChopGetVolume2D(nvecx,nvecy,dist,dx,dy);
        
        thickness = 1*dx;
        % if rr-R < -thickness/2
            % unode(i,j) = 1.0;
        % elseif rr-R > thickness/2
            % unode(i,j) = 0.0;
        % else
            % unode(i,j) = 1.0 - (rr-R+thickness/2) / thickness;
        % end
        unode(i,j) = 0.5*tanh(-(rr-R)/thickness) + 0.5;
    else
        if rr < R
            unode(i,j) = 1.0;
        else
            unode(i,j) = 0.0;
        end
    end
end
end

figure;
contourf(xnode,ynode,unode);
axis equal;

I = 2:nx;
J = 2:ny;
ux = 1.0/(2*dx) .* (unode(I+1,J) - unode(I-1,J));
uy = 1.0/(2*dy) .* (unode(I,J+1) - unode(I,J-1));
uu = sqrt(ux.^2 + uy.^2);

gusum = sum(uu(:)) * dx*dy
scirc = 2*pi*R

usum = sum(unode(:)) * dx*dy
vcirc = pi*R^2






