
function [] = efg2d()

clear;

%%
Lb = 48;
D = 12;

%
young = 30e6;
nu = 0.3;
P = 1000;

% plane stress
Dmat = (young/(1-nu^2)) * [1,nu,0;nu,1,0;0,0,(1-nu)/2];

%
quado = 4;


%% nodes
ndivl = 10;
ndivw = 4;
[x,~,~,numnod] = mesh2(Lb,D,ndivl,ndivw);

% domain of influence
dmax = 3.5;
xspac = Lb / ndivl;
yspac = D / ndivw;
dm(1,1:numnod) = dmax*xspac;
dm(2,1:numnod) = dmax*yspac;

%% quadrature cells
ndivlq = 10;
ndivwq = 4;
[xc,conn,numcell,numq] = mesh2(Lb,D,ndivlq,ndivwq);

%% setup BC
ind1 = 0;
ind2 = 0;
for j = 1:numnod
    % left edge
    if x(1,j) == 0
        ind1 = ind1 + 1;
        nnu(1,ind1) = x(1,j);
        nnu(2,ind1) = x(2,j);
    end
    
    % right edge
    if x(1,j) == Lb
        ind2 = ind2 + 1;
        nt(1,ind2) = x(1,j);
        nt(2,ind2) = x(2,j);
    end
end
lthu = length(nnu);
ltht = length(nt);
ubar = zeros(lthu*2,1);
fext = zeros(numnod*2,1);

% guass points on traction BC
ind = 0;
gauss = pgauss(quado);
for i = 1:(ltht-1)
    ycen = (nt(2,i)+nt(2,i+1))/2;
    jcob = abs(nt(2,i+1)-nt(2,i)) / 2;
    for j = 1:quado
        mark(j) = ycen - gauss(1,j) * jcob;
        ind = ind + 1;
        gst(1,ind) = nt(1,i);
        gst(2,ind) = mark(j);
        gst(3,ind) = gauss(2,j);
        gst(4,ind) = jcob;
    end
end

% gauss points on disp BC
gsu = gst;
gsu(1,:) = 0.0;
qk = zeros(2*lthu,1);


if 1
    figure;
    plot(x(1,:),x(2,:),'.');
    axis equal;
    
    hold on;
    
    % quadrature cells
    for i = 1:numcell
        ii = conn(:,i);
        plot(xc(1,ii),xc(2,ii),'r:');
    end
    
    % gauss points on boundary
    plot(gst(1,:),gst(2,:),'g+', gsu(1,:),gsu(2,:),'mx');
    
    rectangle('Position',[0,-D/2,Lb,D]);
    
    hold off;
end


return
end


function [x,conn,numcell,numq] = mesh2(length,height,ndivl,ndivw)
    
    numcell = ndivl * ndivw;
    numq = (ndivl+1) * (ndivw+1);
    
    dl = length / ndivl;
    dh = height / ndivw;
    
    % coordinates
    for i = 1:ndivl+1
    for j = 1:ndivw+1
        ind = (ndivw+1)*(i-1) + j;
        x(1,ind) = dl * (i-1);
        x(2,ind) = -dh * (j-1) + height/2;
    end
    end
    
    % cell connnectivity
    for j = 1:ndivl
    for i = 1:ndivw
        elemn = (j-1)*ndivw + i;
        nodet(elemn,1) = elemn + (j-1);
        nodet(elemn,2) = nodet(elemn,1) + 1;
        nodet(elemn,3) = nodet(elemn,2) + ndivw + 1;
        nodet(elemn,4) = nodet(elemn,3) - 1;
    end
    end
    
    conn = nodet.';
    
    
    return
end


function [v] = pgauss(k)
    xx1 = sqrt(3/7 + 2/7*sqrt(6/5));
    ww1 = (18-sqrt(30))/36;
    xx2 = sqrt(3/7 - 2/7*sqrt(6/5));
    ww2 = (18+sqrt(30))/36;
    v = [ -xx1, -xx2, xx2, xx1; ww1, ww2, ww2, ww1 ];
    return
end






