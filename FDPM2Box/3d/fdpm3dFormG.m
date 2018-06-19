function [G] = fdpm3dFormG(np,dNx,dNy,dNz)
%fdpm3dFormG
% G is the 9-component full strain-displacement
%

G = zeros(9,np*3);

% xx
G(1,1:3:end) = dNx;

% yy
G(2,2:3:end) = dNy;

% zz
G(3,3:3:end) = dNz;

% xy,yx
G(4,1:3:end) = dNy;
G(5,2:3:end) = dNx;

% yz,zy
G(6,2:3:end) = dNz;
G(7,3:3:end) = dNy;

% zx,xz
G(8,3:3:end) = dNx;
G(9,1:3:end) = dNz;



return
end




