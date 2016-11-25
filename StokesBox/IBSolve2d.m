
function [ax] = IBSolve2d(kx,ky, np,fp,gp,up,vp, nx,ny,dx,dy,nu)

global Emat Hmat;

% distribute to grid
fs = Emat * fp;
gs = Emat * gp;
fs = reshape(fs, nx,ny);
gs = reshape(gs, nx,ny);

% 
[us,vs,ps] = FastSolve2d(kx,ky, fs,gs, nx,ny,dx,dy,nu);

us = reshape(us, nx*ny,1);
vs = reshape(vs, nx*ny,1);
uint = Hmat * us;
vint = Hmat * vs;

ax = [uint; vint];

return
end




