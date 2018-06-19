function [K,G,De] = materialLinElastiMod(E,nu)

K = E / 3 / (1-2*nu);
G = E / 2 / (1+nu);

De = E/(1+nu)/(1-2*nu) * [...
1-nu, nu,   0,      nu; 
nu,   1-nu, 0,      nu;
0,    0,    0.5-nu, 0; 
nu,   nu,   0,      1-nu];

return
end

