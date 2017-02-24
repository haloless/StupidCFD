function [val] = ModSphBesselK(nu,x)
% The modified spherical Bessel function of the 2nd kind
% NOTE the form of the current function is slightly different from common definition,
% the pre-factor is a different constant.
%

% val = sqrt((2/pi)./x) .* besselk(nu+0.5,x);
val = sqrt((pi/2)./x) .* besselk(nu+0.5,x);

return
end

