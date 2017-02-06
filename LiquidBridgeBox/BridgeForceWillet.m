function [F] = BridgeForceWillet(R,H,theta,V,sigma)
% Force by curve fitting in (Willet et al., 2000)
% theta<50, V<0.1
%

% dimensionless volume
Vstar = V / R^3;
% half separation
S = H * 0.5;
%
L = sqrt(V/R);
Splus = S / L;


theta2 = theta^2;
lnV = log(Vstar);
lnV2 = lnV^2;
lnV3 = lnV^3;
lnS = log(Splus);

f1 = (-0.44507 + 0.050832*theta - 1.1466*theta2) + ...
(-0.1119 - 0.000411*theta - 0.1490*theta2) * lnV + ...
(-0.012101 - 0.0036456*theta - 0.01255*theta2) * lnV2 + ...
(-0.0005 - 0.0003505*theta - 0.0029076*theta2) * lnV3;

f2 = (1.9222 - 0.57473*theta - 1.2918*theta2) + ...
(-0.0668 - 0.1201*theta - 0.22574*theta2) * lnV + ...
(-0.0013375 - 0.0068988*theta - 0.01137*theta2) * lnV2;

f3 = (1.268 - 0.01396*theta - 0.23566*theta2) + ...
(0.198 + 0.092*theta - 0.06418*theta2) * lnV + ...
(0.02232 + 0.02238*theta - 0.009853*theta2) * lnV2 + ...
(0.0008585 + 0.001318*theta - 0.00053*theta2) * lnV3;

f4 = (-0.010703 + 0.073776*theta - 0.34742*theta2) + ...
(0.03345 + 0.04543*theta - 0.09056*theta2) * lnV + ...
(0.0018574 + 0.004456*theta - 0.006257*theta2) * lnV2;

Fstar = exp(f1 - f2 * exp(f3*lnS + f4*lnS^2));
F = Fstar * 2*pi*R*sigma;



return
end


