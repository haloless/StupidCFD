
function [ gphix, gphiy ] = CalcGradPhi(phi,nx,ny,dx,dy, bctype,bcval)

bc_neu = 0;
bc_dir = 1;
bc_per = 2;

gphix = zeros(nx+1,ny);
gphiy = zeros(nx,ny+1);

I = 2:nx;
J = 1:ny;
gphix(I,J) = 1.0/dx .* (phi(I,J) - phi(I-1,J));
if bctype(1,1) == bc_neu
    gphix(1,J) = bcval(1,1);
elseif bctype(1,1) == bc_dir
    gphix(1,J) = 1.0/(dx/2) .* (phi(1,J) - bcval(1,1));
elseif bctype(1,1) == bc_per
    gphix(1,J) = 1.0/dx .* (phi(1,J) - phi(nx,J) + (bcval(1,2)-bcval(1,1)));
end
if bctype(1,2) == bc_neu
    gphix(nx+1,J) = bcval(1,2);
elseif bctype(1,2) == bc_dir
    gphix(nx+1,J) = 1.0/(dx/2) .* (bcval(1,1) - phi(nx,J));
elseif bctype(1,2) == bc_per
    gphix(nx+1,J) = 1.0/dx .* (phi(1,J) - phi(nx,J) + (bcval(1,2)-bcval(1,1)));
end

I = 1:nx;
J = 2:ny;
gphiy(I,J) = 1.0/dy .* (phi(I,J)-phi(I,J-1));
if bctype(2,1) == bc_neu
    gphiy(I,1) = bcval(2,1);
elseif bctype(2,1) == bc_dir
    gphiy(I,1) = 1.0/(dy/2) .* (phi(I,1) - bcval(2,1));
elseif bctype(2,1) == bc_per
    gphiy(I,1) = 1.0/dy .* (phi(I,1) - phi(I,ny) + (bcval(2,2)-bcval(2,1)));
end
if bctype(2,2) == bc_neu
    gphiy(I,ny+1) = bcval(2,2);
elseif bctype(2,2) == bc_dir
    gphiy(I,ny+1) = 1.0/(dy/2) .* (bcval(2,2) - phi(I,ny));
elseif bctype(2,2) == bc_per
    gphiy(I,ny+1) = 1.0/dy .* (phi(I,1) - phi(I,ny) + (bcval(2,2)-bcval(2,1)));
end


return
end

