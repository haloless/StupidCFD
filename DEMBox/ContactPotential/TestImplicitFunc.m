
clear all;


fun1 = @(x,y) (x/4).^8 + (y/3).^8 - 1;

fun2 = @(x,y) ((x/4).^8 + (y/3).^8).^(1/8) - 1;


dx = 0.1;
dy = 0.1;
[xgrid,ygrid] = ndgrid(-6:dx:6,-6:dy:6);
nx = size(xgrid,1);
ny = size(xgrid,2);

phi1 = fun1(xgrid,ygrid);
phi2 = fun2(xgrid,ygrid);

if 0
figure;
contour(xgrid,ygrid,phi1,'ShowText','on');
hold on;
contour(xgrid,ygrid,phi1,[0,0],'k');
hold off;
end

% bandwidth = 1.0;
% bandwidth = 3.0;
% bandwidth = 4.5;
bandwidth = 5.0;
% bandwidth = 6.0;
% bandwidth = 8.0;

phi1reg = ImplicitFuncReinit(phi1,nx,ny,dx,dy,bandwidth);

figure;
contourf(xgrid,ygrid,phi1reg,'ShowText','on');
hold on;
contour(xgrid,ygrid,phi1,[0,0],'r');
contour(xgrid,ygrid,phi1reg,[0,0],'k');
hold off;
axis equal;

gphix = zeros(nx,ny);
gphiy = zeros(nx,ny);
for i = 1:nx
for j = 1:ny
	if i == 1
		gphix(i,j) = (phi1reg(i+1,j)-phi1reg(i,j)) / dx;
	elseif i == nx
		gphix(i,j) = (phi1reg(i,j)-phi1reg(i-1,j)) / dx;
	else
		gphix(i,j) = (phi1reg(i+1,j)-phi1reg(i-1,j)) / (dx*2);
	end
	
	if j == 1
		gphiy(i,j) = (phi1reg(i,j+1)-phi1reg(i,j)) / dy;
	elseif j == ny
		gphiy(i,j) = (phi1reg(i,j)-phi1reg(i,j-1)) / dy;
	else
		gphiy(i,j) = (phi1reg(i,j+1)-phi1reg(i,j-1)) / (dy*2);
	end
end
end

gphi = sqrt(gphix.^2 + gphiy.^2);

figure;
contourf(xgrid,ygrid,gphi);
colorbar;
% hold on;
% contour(xgrid,ygrid,phi1,[0,0],'r');
% contour(xgrid,ygrid,phi1reg,[0,0],'k');
% hold off;
axis equal;



% figure;
% contour(xgrid,ygrid,phi2,'ShowText','on');
% hold on;
% contour(xgrid,ygrid,phi2,[0,0],'k');
% hold off;


