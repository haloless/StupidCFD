
clc;
clear all;

xoffset = 3.0;
yoffset = 4.0;
zoffset = 5.0;

R = 1.0;
V = 4/3*pi*R^3;
Isphere = 2/5*V*R^2;
nsub = 80;
hsub = R*2 / nsub;
vsub = hsub^3;

gridx = linspace(-R+hsub/2,R-hsub/2,nsub);
gridy = linspace(-R+hsub/2,R-hsub/2,nsub);
gridz = linspace(-R+hsub/2,R-hsub/2,nsub);

[gridx,gridy,gridz] = ndgrid(gridx,gridy,gridz);
gridx = gridx + xoffset;
gridy = gridy + yoffset;
gridz = gridz + zoffset;

mask = ((gridx-xoffset).^2 + (gridy-yoffset).^2 + (gridz-zoffset).^2 <= R^2);

Ixy = -sum(gridx(mask).*gridy(mask)) * vsub;
Iyz = -sum(gridy(mask).*gridz(mask)) * vsub;
Izx = -sum(gridz(mask).*gridx(mask)) * vsub;

Ixy0 = -xoffset * yoffset * V;
Iyz0 = -yoffset * zoffset * V;
Izx0 = -zoffset * xoffset * V;

Ixx = sum(gridy(mask).^2 + gridz(mask).^2) * vsub;
Iyy = sum(gridz(mask).^2 + gridx(mask).^2) * vsub;
Izz = sum(gridx(mask).^2 + gridy(mask).^2) * vsub;

Ixx0 = (yoffset^2+zoffset^2)*V + Isphere;
Iyy0 = (zoffset^2+xoffset^2)*V + Isphere;
Izz0 = (xoffset^2+yoffset^2)*V + Isphere;

J0 = [Ixx0,Ixy0,Izx0; Ixy0,Iyy0,Iyz0; Izx0,Iyz0,Izz0];

Jder = DeriveShiftedMOI(0.0,0.0,0.0, 1, xoffset,yoffset,zoffset, [1.0*V], Isphere*eye(3));

