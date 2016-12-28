function [muf] = PartWallMobUF(theta,ah,ewall)

[f1,f3,f5,g1,g3,g5] = InterpUF(theta);

ah3 = ah^3;
ah5 = ah^5;

cf = 1 - f1*ah + f3*ah3 - f5*ah5;
cg = 1 - g1*ah + g3*ah3 - g5*ah5;

ww = ewall * ewall';
muf = (eye(3)-ww).*cf + ww.*cg;

return
end


