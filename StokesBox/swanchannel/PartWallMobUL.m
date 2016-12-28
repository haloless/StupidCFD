function [mul] = PartWallMobUL(theta,ah,ewall)

[f2,f4] = InterpUL(theta);

ah2 = ah^2;
ah4 = ah^4;

cf = f2*ah2 + f4*ah4;

cf = 1 - f1*ah + f3*ah3 - f5*ah5;
cg = 1 - g1*ah + g3*ah3 - g5*ah5;

ww = ewall * ewall';
muf = (eye(3)-ww).*cf + ww.*cg;

return
end


