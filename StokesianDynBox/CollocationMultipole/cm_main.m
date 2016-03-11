
% clc;
clear all;

% collocation points
Np = 12
% Np = 48
% center-to-center sphere separation
R = 4.0

rfunc = {};

% set collocation points
ang1 = zeros(Np,1);
cos1 = zeros(Np,1);
for i = 1:Np
    theta = (i-1) * pi/(Np-1);
    ang1(i) = theta;
    cos1(i) = cos(theta);
end

%
% X11A, X12A, X11G, X12G
%
rhs = zeros(Np*3,1);
for i = 1:Np
    i2 = i + Np;
    i3 = i2 + Np;
    rhs(i) = 1.0;
    rhs(i2) = 0.0;
    rhs(i3) = 0.0;
end
% R,nstar,m,sgn,np,rhs
[ amn, bmn, cmn ] = cm_stokes(R,1,0,-1,Np,rhs);
t1 = 2.0/3.0 * amn(1);
t3 = -0.5 * amn(2);

for i = 1:Np
    i2 = i + Np;
    i3 = i2 + Np;
    rhs(i) = 1.0;
    rhs(i2) = 0.0;
    rhs(i3) = 0.0;
end
[ amn, bmn, cmn ] = cm_stokes(R,1,0,1,Np,rhs);
t2 = 2.0/3.0 * amn(1);
t4 = -0.5 * amn(2);

rfunc.x11a = 0.5 * (t1 + t2);
rfunc.x12a = 0.5 * (t1 - t2);
rfunc.x11g = 0.5 * (t3 + t4);
rfunc.x12g = 0.5 * (t3 - t4);


%
% Y11A, Y12A, Y11B, Y12B, Y11G, Y12G
%
rhs = zeros(Np*3,1);
for i = 1:Np
    i2 = i + Np;
    i3 = i2 + Np;
    rhs(i) = 0.0;
    rhs(i2) = 0.0;
    rhs(i3) = 1.0;
end

[ amn, bmn, cmn ] = cm_stokes(R,1,1,-1,Np,rhs);
t1= 2.0/3.0 * amn(1);
t3 = 2.0 * cmn(1);
t5 = -0.5 * amn(2);

[ amn, bmn, cmn ] = cm_stokes(R,1,1,1,Np,rhs);
t2 = 2.0/3.0 * amn(1);
t4 = 2.0 * cmn(1);
t6 = -0.5 * amn(2);

rfunc.y11a = 0.5 * (t1 + t2);
rfunc.y12a = 0.5 * (t2 - t1);
rfunc.y11b = -0.5 * (t3 + t4);
rfunc.y12b = 0.5 * (t3 - t4);
rfunc.y11g = 0.5 * (t5 + t6);
rfunc.y12g = 0.5 * (t6 - t5);

%
% X11C, X12C
%
rhs = zeros(Np,1);
for i = 1:Np
    rhs(i) = -1.0;
end
cmn = cm_stokesd(R,-1,Np,rhs);
t1 = cmn(1);

for i = 1:Np
    rhs(i) = -1.0;
end
cmn = cm_stokesd(R,1,Np,rhs);
t2 = cmn(1);

rfunc.x11c = 0.5 * (t1 + t2);
rfunc.x12c = -0.5 * (t1 - t2);


%
% Y11C, Y12C, Y11H, Y12H
%
rhs = zeros(Np*3,1);
for i = 1:Np
    i2 = i + Np;
    i3 = i2 + Np;
    rhs(i) = 1.0;
    rhs(i2) = 0.0;
    rhs(i3) = -cos1(i);
end

[ amn, bmn, cmn ] = cm_stokes(R,1,1,-1,Np,rhs);
t1 = 2.0/3.0 * amn(1);
t3 = cmn(1);
t5 = -0.25 * amn(2);

[ amn, bmn, cmn ] = cm_stokes(R,1,1,1,Np,rhs);
t2 = 2.0/3.0 * amn(1);
t4 = cmn(1);
t6 = -0.25 * amn(2);

rfunc.y11c = 0.5 * (t3 + t4);
rfunc.y12c = 0.5 * (t3 - t4);
rfunc.y11h = 0.5 * (t5 + t6);
rfunc.y12h = 0.5 * (t5 - t6);


%
% X11+X12M, Y11+Y12M, Z11+Z12M
%
rhs = zeros(Np*3,1);
for i = 1:Np
    i2 = i + Np;
    i3 = i2 + Np;
    rhs(i) = 2.0 * cos1(i);
    rhs(i2) = -1.0;
    rhs(i3) = 0.0;
end

[ amn, bmn, cmn ] = cm_stokes(R,1,0,1,Np,rhs);
rfunc.x1112m = 0.1 * amn(2);

rhs = zeros(Np*3,1);
for i = 1:Np
    i2 = i + Np;
    i3 = i2 + Np;
    rhs(i) = 3.0;
    rhs(i2) = 0.0;
    rhs(i3) = 3.0 * cos1(i);
end
[ amn, bmn, cmn ] = cm_stokes(R,1,1,-1,Np,rhs);
rfunc.y1112m = 0.1 * amn(2);

rhs = zeros(Np*3,1);
for i = 1:Np
    i2 = i + Np;
    i3 = i2 + Np;
    rhs(i) = 0.0;
    rhs(i2) = 0.0;
    rhs(i3) = 6.0;
end
[ amn, bmn, cmn ] = cm_stokes(R,2,2,1,Np,rhs);
rfunc.z1112m = 0.1 * amn(1);


rfunc
