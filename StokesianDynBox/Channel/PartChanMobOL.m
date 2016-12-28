function [mol] = PartChanMobOL(theta,ah,ewall)

[f3,g3] = InterpOL(theta);

ah3 = ah^3;

cf = 0.75 - f3*ah3;
cg = 0.75 - g3*ah3;

ww = ewall * ewall';
mol = (eye(3)-ww).*cf + ww.*cg;

return
end


