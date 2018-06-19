
clear;

ndim = 2;

Lx = 1;
Ly = 1;

% nx = 10;
% nx = 20;
% nx = 40;
% nx = 80;
nx = 100;

dh = Lx/nx;
ny = round(Ly/dh);

re = dh * 2.1;
% re = dh * 3.1;

wfunc = @(r) (re./r - 1);

% calculate n0
nrange = round(re/dh+1);
nzero = 0;
for i = -nrange:nrange
for j = -nrange:nrange
    if i~=0 || j~=0
        rr = sqrt((i*dh)^2 + (j*dh)^2);
        if rr < re
            nzero = nzero + wfunc(rr);
        end
    end
end
end


% generate particle
partnum = 0;
xpart = [];
ypart = [];
volpart = [];
surfflag = [];
for j = 1:ny
for i = 1:nx
    xx = (i-0.5)*dh;
    yy = (j-0.5)*dh;
    
    if j <= i
        
        partnum = partnum + 1;
        xpart(end+1,1) = xx;
        ypart(end+1,1) = yy;
        volpart(end+1,1) = dh^2;
        surfflag(end+1,1) = logical(j==i || j==i-1);
        
    end
end
end
surfflag = logical(surfflag);
surflist = find(surfflag)';
nsurf = numel(surflist);





npart = [];
for ipart = 1:partnum
    isurf = ipart;
    rx = xpart - xpart(ipart);
    ry = ypart - ypart(ipart);
    rr = sqrt(rx.^2 + ry.^2);
    reflag = (rr<re); reflag(ipart) = 0;
    relist = find(reflag);
    rewgt = wfunc(rr(relist));
    
    npart(isurf,1) = sum(rewgt);
    
    conn(isurf).neighnum = numel(relist);
    conn(isurf).neighlist = relist';
    conn(isurf).neighwgt = rewgt;
    conn(isurf).neighrx = rx(relist);
    conn(isurf).neighry = ry(relist);
    conn(isurf).neighrr = rr(relist);
end


% surfflag = (npart<0.95*nzero);
% surfflag(xpart.^2+ypart.^2<(Ra+Rb)/2) = 0;
% surflist = find(surfflag)';
% nsurf = numel(surflist);

if 1
    figure;
    plot(xpart,ypart,'o', xpart(surfflag),ypart(surfflag),'x');
    axis equal;
    % axis([0 Lx 0 Ly]);
    hold on;
    rectangle('Position',[0 0 Lx Ly]);
    hold off;
    % return
end
if 0
    figure;
    for isurf = 1:nsurf
        plot(xpart,ypart,'o');
        hold on;
        ipart = surflist(isurf);
        plot(xpart(ipart),ypart(ipart),'x', xpart(conn(isurf).neighlist),ypart(conn(isurf).neighlist),'+');
        hold off;
        axis equal;
        pause;
    end
end



% qxfunc = @(x,y) -x ./ sqrt(x.^2+y.^2);
% qyfunc = @(x,y) -y ./ sqrt(x.^2+y.^2);
% qdivfunc = @(x,y) -1 ./ sqrt(x.^2+y.^2);
% qxfunc = @(x,y) zeros(size(x));
qxfunc = @(x,y) 2*ones(size(x));
qyfunc = @(x,y) -ones(size(y));
qdivfunc = @(x,y) zeros(size(x));

spart = [];
if 0
    for isurf = 1:nsurf
        ipart = surflist(isurf);
        ineigh = conn(ipart).neighlist;
        ineighrx = conn(ipart).neighrx;
        ineighry = conn(ipart).neighry;
        ineighr2 = conn(ipart).neighrr.^2;
        ineighw = conn(ipart).neighwgt;
        
        qxi = qxfunc(xpart(ipart),ypart(ipart));
        qyi = qyfunc(xpart(ipart),ypart(ipart));
        
        qxs = qxfunc(xpart(ineigh),ypart(ineigh));
        qys = qyfunc(xpart(ineigh),ypart(ineigh));
        
        qqs = (qxs.*ineighrx+qys.*ineighry) ./ ineighr2 .* ineighw;
        
        divinner = ndim/nzero * sum(qqs(~surfflag(ineigh)));
        
        divqext = qdivfunc(xpart(ipart),ypart(ipart));
        
        spart(isurf,1) = -divqext + divinner;
        % spart(isurf,1) = -divqext;
    end
end
if 0
    for isurf = 1:nsurf
        ipart = surflist(isurf);
        ineigh = conn(ipart).neighlist;
        ineighrx = conn(ipart).neighrx;
        ineighry = conn(ipart).neighry;
        ineighr2 = conn(ipart).neighrr.^2;
        ineighw = conn(ipart).neighwgt;
        
        qxi = qxfunc(xpart(ipart),ypart(ipart));
        qyi = qyfunc(xpart(ipart),ypart(ipart));
        
        qxs = qxfunc(xpart(ineigh),ypart(ineigh));
        qys = qyfunc(xpart(ineigh),ypart(ineigh));
        qxs(~surfflag(ineigh)) = 0;
        qys(~surfflag(ineigh)) = 0;
        
        qqs = ((qxs-qxi).*ineighrx + (qys-qyi).*ineighry) ./ ineighr2 .* ineighw;
        
        spart(isurf,1) = -ndim/nzero * sum(qqs);
    end
end
if 1
    % mid-particle
    for isurf = 1:nsurf
        ipart = surflist(isurf);
        ineigh = conn(ipart).neighlist;
        ineighrx = conn(ipart).neighrx;
        ineighry = conn(ipart).neighry;
        ineighr2 = conn(ipart).neighrr.^2;
        ineighw = conn(ipart).neighwgt;
        
        qxi = qxfunc(xpart(ipart),ypart(ipart));
        qyi = qyfunc(xpart(ipart),ypart(ipart));
        
        qxs = qxfunc(xpart(ineigh),ypart(ineigh));
        qys = qyfunc(xpart(ineigh),ypart(ineigh));
        % qxs = (qxs + qxi) / 2;
        % qys = (qys + qxi) / 2;
        % qxs = qxs - qxi;
        % qys = qys - qyi;
        % qxs = qxs + qxi;
        % qys = qys + qyi;
        
        qqs = (qxs.*ineighrx+qys.*ineighry) ./ ineighr2 .* ineighw;
        
        divinner = 2*ndim/nzero * sum(qqs(:));
        
        divqext = qdivfunc(xpart(ipart),ypart(ipart));
        
        spart(isurf,1) = -divqext + divinner;
        % spart(isurf,1) = -divqext;
    end
end
svolpart = spart .* volpart(surfflag);
sum(svolpart)

if 0
    ipart = floor(partnum/2);
    assert(~surfflag(ipart));
    
    ineigh = conn(ipart).neighlist;
    ineighrx = conn(ipart).neighrx;
    ineighry = conn(ipart).neighry;
    ineighr2 = conn(ipart).neighrr.^2;
    ineighw = conn(ipart).neighwgt;
    
    qxi = qxfunc(xpart(ipart),ypart(ipart));
    qyi = qyfunc(xpart(ipart),ypart(ipart));
    
    qxs = qxfunc(xpart(ineigh),ypart(ineigh));
    qys = qyfunc(xpart(ineigh),ypart(ineigh));
    
    divq = ((qxs-qxi).*ineighrx + (qys-qyi).*ineighry) ./ ineighr2 .* ineighw;
    divq = ndim/nzero * sum(divq);
    
    [divq, qdivfunc(xpart(ipart),ypart(ipart))]
end













