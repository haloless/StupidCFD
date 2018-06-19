function [wp,xsi,eta,zet,dNr] = der_shape_func(ngp)
% Gauss points and Shape function derivatives


sqrt3 = sqrt(3);

xsi = [-1; -1;  1;  1; -1; -1;  1;  1] ./ sqrt3;
eta = [-1; -1; -1; -1;  1;  1;  1;  1] ./ sqrt3;
zet = [-1;  1   1; -1; -1;  1;  1; -1] ./ sqrt3;

wp = ones(8,1);

r2 = ngp * 3;

dNr(1:3:r2,1) = -1/8 .* (1-eta) .* (1-zet);
dNr(1:3:r2,2) = -1/8 .* (1-eta) .* (1+zet);
dNr(1:3:r2,3) =  1/8 .* (1-eta) .* (1+zet);
dNr(1:3:r2,4) =  1/8 .* (1-eta) .* (1-zet);
dNr(1:3:r2,5) = -1/8 .* (1+eta) .* (1-zet);
dNr(1:3:r2,6) = -1/8 .* (1+eta) .* (1+zet);
dNr(1:3:r2,7) =  1/8 .* (1+eta) .* (1+zet);
dNr(1:3:r2,8) =  1/8 .* (1+eta) .* (1-zet);

dNr(2:3:r2+1,1) = -1/8 .* (1-xsi) .* (1-zet);
dNr(2:3:r2+1,2) = -1/8 .* (1-xsi) .* (1+zet);
dNr(2:3:r2+1,3) = -1/8 .* (1+xsi) .* (1+zet);
dNr(2:3:r2+1,4) = -1/8 .* (1+xsi) .* (1-zet);
dNr(2:3:r2+1,5) =  1/8 .* (1-xsi) .* (1-zet);
dNr(2:3:r2+1,6) =  1/8 .* (1-xsi) .* (1+zet);
dNr(2:3:r2+1,7) =  1/8 .* (1+xsi) .* (1+zet);
dNr(2:3:r2+1,8) =  1/8 .* (1+xsi) .* (1-zet);

dNr(3:3:r2+2,1) = -1/8 .* (1-xsi) .* (1-eta);
dNr(3:3:r2+2,2) =  1/8 .* (1-xsi) .* (1-eta);
dNr(3:3:r2+2,3) =  1/8 .* (1+xsi) .* (1-eta);
dNr(3:3:r2+2,4) = -1/8 .* (1+xsi) .* (1-eta);
dNr(3:3:r2+2,5) = -1/8 .* (1-xsi) .* (1+eta);
dNr(3:3:r2+2,6) =  1/8 .* (1-xsi) .* (1+eta);
dNr(3:3:r2+2,7) =  1/8 .* (1+xsi) .* (1+eta);
dNr(3:3:r2+2,8) = -1/8 .* (1+xsi) .* (1+eta);



return
end


