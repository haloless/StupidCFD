function [mul] = PartChanMobUL(theta,ah,ewall)

CommonGlobals;

[f2,f4] = InterpUL(theta);

ah2 = ah^2;
ah4 = ah^4;

cf = f2*ah2 + f4*ah4;

mul = zeros(3,3);
for i = 1:3
for j = 1:3
	for k = 1:3
		mul(i,j) = mul(i,j) + Epsil(i,j,k)*ewall(k);
	end
end
end

mul = mul.*cf;

return
end


