
clear all;

EKConst;
EKInitConst;

H = 1.0;
nx = 128;
% H = 16.0;
% nx = 80;
% H = 128.0;
% nx = 160;

xlo = -H / 2;
xhi = H / 2;
dx = (xhi-xlo) / nx;
dx2 = dx * dx;
xcell = linspace(xlo+dx/2,xhi-dx/2,nx)';
xnode = linspace(xlo,xhi,nx+1)';

% boundary potential
psi0 = 1.0;
psi_xlo = psi0;
psi_xhi = psi0;
% Debye length
kappa = 1.0 / (H*0.01);
kappa2 = kappa^2;


% build 1D Laplacian
Lap = sparse(nx,nx);
rbc = zeros(nx,1);
for i = 1:nx
    if i > 1
        Lap(i,i-1) = 1.0/dx2;
        Lap(i,i) = Lap(i,i) - 1.0/dx2;
    else
        Lap(i,i) = Lap(i,i) - 2.0/dx2;
        rbc(i) = rbc(i) - 2.0*psi_xlo/dx2;
    end
    if i < nx
        Lap(i,i+1) = 1.0/dx2;
        Lap(i,i) = Lap(i,i) - 1.0/dx2;
    else
        Lap(i,i) = Lap(i,i) - 2.0/dx2;
        rbc(i) = rbc(i) - 2.0*psi_xhi/dx2;
    end
end

% electric potential
psi = zeros(nx,1);

if (1)
    % linearized Poisson-Boltzmann
    % as initial guess
    D = kappa2 .* speye(nx);
    L = Lap - D;
    rhs = rbc;
    psi = L \ rhs;
end

if (1)
    conv = 0;
    ncycle = 100;
    omega = 1.0;
    tol_rel = 1.0e-14;
    resid = Lap*psi - rbc - kappa2.*sinh(psi);
    rnorm0 = norm(resid);
    tol = tol_rel * rnorm0;
    
    for cycle = 1:ncycle
        psi_old = psi;
        
        diag = kappa2 .* cosh(psi_old);
        D = spdiags(diag,0,nx,nx);
        
        rhs = kappa2.*sinh(psi_old) - diag.*psi_old;
        rhs = rhs + rbc;
        
        L = Lap - D;
        
        psi = L \ rhs;
        
        % under-relaxation
        psi = omega.*psi + (1.0-omega).*psi_old;
        
        % calculate non-linear residual
        resid = Lap*psi - rbc - kappa2.*sinh(psi);
        rnorm = norm(resid);
        
        % disp(['cycle=',int2str(cycle),';rnorm=',num2str(rnorm)]);
        disp(['cycle=',int2str(cycle),';r/r0=',num2str(rnorm/rnorm0)]);
        
        if (rnorm < tol) 
            conv = 1;
            break;
        end
    end
    
    if (~conv)
        disp(['Failed']);
    end
end


if (1)
    % linearized analytical solution
    psi0 = psi_xlo;
    psi_ana = psi0 .* cosh(kappa.*xnode) ./ cosh(0.5*kappa*H);
    
    figure;
    plot(xnode,psi_ana,'-', xcell,psi,'x');
    axis([xlo,xhi,0,max(psi_xlo,psi_xhi)]);
    legend('linear','numerical');
end

