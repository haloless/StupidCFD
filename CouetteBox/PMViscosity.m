


function [ nucell nunode ] = PMVelocityLapOp (eb_vof,nx,ny)
%
EBGlobals;

% compute (penalty) viscosity from harmonic average

% cell centered
nucell = zeros(nx+2,ny+2);
nucell = nu*nu_s ./ (eb_vof.*nu + (1.0-eb_vof).*nu_s);

% node centered
nunode = zeros(nx+1,ny+1);
% compute node fraction first
I = 1:nx+1;
J = 1:ny+1;
vofnode = 0.25 * (eb_vof(I,J) + eb_vof(I+1,J) + eb_vof(I,J+1) + eb_vof(I+1,J+1));
%
nunode = nu*nu_s ./ (vofnode.*nu + (1.0-vofnode).*nu_s);


return
end



