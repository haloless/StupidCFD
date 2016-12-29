
function [ Pmat ] = InductionMat(a,h,H)

theta = h / H;

% Pmat = zeros(2,2);

euler_const = 0.5772156649;
zeta3 = 1.202056903159594;

p11 = 1 + 0.5*(psi(0,theta) + psi(0,1-theta) + 2*euler_const) * (a/H);
p12 = 0.25 * (psi(1,theta) - psi(1,1-theta)) * (a/H)^2;
p22 = 1 + 0.125 * (psi(2,theta) + psi(2,1-theta) - 4*zeta3) * (a/H)^3;


Pmat = [ p11,p12; p12,p22 ];


return
end


