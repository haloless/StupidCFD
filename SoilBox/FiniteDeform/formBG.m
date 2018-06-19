function [B,G] = formBG(dNx,nen)
% B is the 6-component standard strain-displacement
% G is the 9-component full strain-displacement


ix = 1:3:nen*3;
iy = 2:3:nen*3;
iz = 3:3:nen*3;

B = zeros(6,nen*3);
% xx
B(1,1:3:end) = dNx(1,:);
% yy
B(2,2:3:end) = dNx(2,:);
% zz
B(3,3:3:end) = dNx(3,:);
% xy
B(4,1:3:end) = dNx(2,:);
B(4,2:3:end) = dNx(1,:);
% yz
B(5,2:3:end) = dNx(3,:);
B(5,3:3:end) = dNx(2,:);
% zx
B(6,1:3:end) = dNx(3,:);
B(6,3:3:end) = dNx(1,:);

G = zeros(9,nen*3);
% xx,yy,zz
G(1:3,:) = B(1:3,:);
% xy,yx
G(4,1:3:end) = dNx(2,:);
G(5,2:3:end) = dNx(1,:);
% yz,zy
G(6,2:3:end) = dNx(3,:);
G(7,3:3:end) = dNx(2,:);
% zx,xz
G(8,3:3:end) = dNx(1,:);
G(9,1:3:end) = dNx(3,:);




return
end

