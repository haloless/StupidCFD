
function LBM_PoiseuilleMain

refine = 1;
Lx = 64*refine;
Ly = 32*refine;
nx = Lx;
ny = Ly;
Nx = nx + 2;
Ny = ny + 2;

hx = Lx / nx;
hy = Ly / ny;
hh = min([hx hy]);
qc = 1;
dt = hh / qc;
cs = 1/sqrt(3) * qc;
cs2 = 1/3 * qc^2;

u = zeros(Nx,Ny);
v = zeros(Nx,Ny);
feq = zeros(Nx,Ny,9);
f = zeros(Nx,Ny,9);




return
end % end of main function




