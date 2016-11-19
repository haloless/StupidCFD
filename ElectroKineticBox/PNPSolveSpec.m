
function [ sol,flag,relres,iter ] = PNPSolveSpec(c,z, gphix,gphiy, sdf,tag, nx,ny,dx,dy, bctype)
% solve species concentration
%

csum0 = sum(c(tag==1));


[flagx,flagy] = NPFlag(sdf,tag,nx,ny,bctype);

afunc = @(x) NPApply(x,z, gphix,gphiy, flagx,flagy, sdf,tag, nx,ny,dx,dy,bctype);

sol = reshape(c, nx*ny,1);
rhs = zeros(nx,ny);
if (0)
% far-field
rhs(:,ny) = 1.0;
end
rhs = reshape(rhs,nx*ny,1);
if (1)
% conservation
rhs(1) = csum0;
end

[sol,flag,relres,iter] = bicgstab(afunc, rhs, 1e-8, 2000, [],[],sol);

% cnew = reshape(sol,nx,ny);

% enforce conservation
% cnt = sum(tag==1);
% csum = sum(cnew(tag==1));
% cerr = csum0 - csum;
% cnew = cnew + cerr/cnt;

return
end



function [ flagx,flagy ] = NPFlag(sdf,tag,nx,ny,bctype)
% create flag for flux
flagx = zeros(nx+1,ny);
flagy = zeros(nx,ny+1);

I = 2:nx; J = 1:ny;
flagx(I,J) = (tag(I-1,J)==1) & (tag(I,J)==1);
if bctype(1,1) == 2
    flagx(1,J) = 1;
end
if bctype(1,2) == 2
    flagx(nx+1,J) = 1;
end

I = 1:nx; J = 2:ny;
flagy(I,J) = (tag(I,J-1)==1) & (tag(I,J)==1);
if bctype(2,1) == 2
    flagy(I,1) = 1;
end
if bctype(2,2) == 2
    flagy(I,ny+1) = 1;
end

return
end


function [ fx,fy ] = NPFlux(c,z, gphix,gphiy, flagx,flagy, nx,ny,dx,dy, bctype)
% flux
fx = zeros(nx+1,ny);
fy = zeros(nx,ny+1);

I = 2:nx; J = 1:ny;
fx(I,J) = 1.0/dx.*(c(I,J)-c(I-1,J)) + 0.5*z.*(c(I,J)+c(I-1,J)).*gphix(I,J);
if bctype(1,1) == 2
    fx(1,J) = 1.0/dx.*(c(1,J)-c(nx,J)) + 0.5*z.*(c(1,J)+c(nx,J)).*gphix(1,J);
end
if bctype(1,2) == 2
    fx(nx+1,J) = 1.0/dx.*(c(1,J)-c(nx,J)) + 0.5*z.*(c(1,J)+c(nx,J)).*gphix(nx+1,J);
end

I = 1:nx; J = 2:ny;
fy(I,J) = 1.0/dy.*(c(I,J)-c(I,J-1)) + 0.5*z.*(c(I,J)+c(I,J-1)).*gphiy(I,J);
if bctype(2,1) == 2
    fy(I,1) = 1.0/dy.*(c(I,1)-c(I,ny)) + 0.5*z.*(c(I,1)+c(I,ny)).*gphiy(I,1);
end
if bctype(2,2) == 2
    fy(I,ny+1) = 1.0/dy.*(c(I,1)-c(I,ny)) + 0.5*z.*(c(I,1)+c(I,ny)).*gphiy(I,ny+1);
end

fx = fx .* flagx;
fy = fy .* flagy;
return
end

function [ q ] = NPApply(x,z, gphix,gphiy, flagx,flagy, sdf,tag, nx,ny,dx,dy,bctype)

c = reshape(x,nx,ny);

[fx,fy] = NPFlux(c,z, gphix,gphiy, flagx,flagy, nx,ny,dx,dy,bctype);

I = 1:nx; J = 1:ny;
q = 1.0/dx.*(fx(I+1,J)-fx(I,J)) + 1.0/dy.*(fy(I,J+1)-fy(I,J));

% 
q(:,ny) = c(:,ny);

q = reshape(q,nx*ny,1);

if 0
    csum = sum(c(tag==1));
    q(1) = csum;
end

return
end

