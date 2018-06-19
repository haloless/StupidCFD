function [p,s] = voigt3dPressShear(vsig)
%voigt3dPressShear

bm1 = [1;1;1;0;0;0];

% volumetric pressure
p = mean(vsig(1:3));

% deviatoric shear
s = vsig - p.*bm1;


return
end

