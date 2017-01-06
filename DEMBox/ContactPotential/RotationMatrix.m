function [rotmat] = RotationMatrix(rotang)
crot = cos(rotang);
srot = sin(rotang);
rotmat = [crot,-srot;srot,crot];

return
end


