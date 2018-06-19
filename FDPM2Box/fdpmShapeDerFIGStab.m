function [NxStab,NyStab] = fdpmShapeDerFIGStab(Nx,Ny,Nxx,Nxy,Nyy, hx,hy)
%fdpmShapeDerFIGStab
% Apply Finite-Increment-Gradient Stabilization to derivative shape functions
% Inputs
% hx,hy: controls the increment size
%

NxStab = Nx + hx*Nxx + hy*Nxy;
NyStab = Ny + hy*Nyy + hx*Nxy;
% NxStab = Nx - hx*Nxx - hy*Nxy;
% NyStab = Ny - hy*Nyy - hx*Nxy;





return
end

