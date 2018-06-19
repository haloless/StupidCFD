function [lambda0,mu0] = materialLameConst(E0,nu0)
%material_LameConst
% Calculate Lame's constants (lambda,mu) from Young E and Poisson nu



lambda0 = E0*nu0 / (1+nu0) / (1-2*nu0);
mu0 = E0 / 2 / (1+nu0);

return
end

