function [B] = fdpm3dFormB(np,dNx,dNy,dNz)
%fdpm3dFormB
% B is the 6-component standard strain-displacement
%

B = zeros(6,np*3);

% xx
B(1,1:3:end) = dNx;
% yy
B(2,2:3:end) = dNy;
% zz
B(3,3:3:end) = dNz;
% xy
B(4,1:3:end) = dNy;
B(4,2:3:end) = dNx;
% yz
B(5,2:3:end) = dNz;
B(5,3:3:end) = dNy;
% zx
B(6,1:3:end) = dNz;
B(6,3:3:end) = dNx;



return
end




