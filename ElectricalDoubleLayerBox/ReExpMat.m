function [Tmat,Rmat,Smat] = ReExpMat(psrc,pdst, kappa,nmax)
% Wrapper function for re-expansion operation: src -> dst

ImUnit = 1i;

% reexpansion direction
pab = pdst - psrc;

% create original->coaxial rotation
ReExpSimpleRot;

% distance
[r,theta,phi] = sh_cart2sph(pab);
kr = kappa * r;

if 1
    theta = thetaprime;
    phi = phiprime;
end

% cutoff
nmax1 = nmax + 1;
npole = nmax1^2;
% disp(['nmax=',int2str(nmax), '; npole=',int2str(npole)]);

% rotation matrix R
ReExpCalcR;

% scale matrix S
ReExpCalcS;

% combine total operations
Tmat = Rmat * Smat * Rmat';

return
end


