
clear;

I = voigt_eye();

% set material 
matl = struct();
matl.nu = 0.3;
matl.E = 25000;
matl.alpha0 = 0.7;
matl.beta0 = 0.7;
matl.a = 0.25;
matl.k = 0.1;

% derive others
% bulk modulus
matl.K = matl.E / (3*(1-2*matl.nu));
% Lame const
matl.mu = matl.E / (2*(1+matl.nu));
matl.lambda = matl.E*matl.nu/((1+matl.nu)*(1-2*matl.nu));
% elastic tensor 6x6
matl.Ce = 2*matl.mu*eye(6) + matl.lambda*(I*I');
% elastic matrix in principal axis
matl.ae = (matl.K + 4*matl.mu/3)*eye(3) + (mal.K-2*matl.mu/3)*(~eye(3));






% initial stress, hydrostatic compression
sigma = -50.0 * I;






