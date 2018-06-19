
clear;

ndim = 2;

Ra = 1;
Rb = 2;

nrad = 10;
% nrad = 20;
% nrad = 40;
% nrad = 80;
% nrad = 100;
% dh = (Rb-Ra)/nrad;
dh = (Rb)/nrad;

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
% surfflag = [];
npart = [];

if 0
    % circular align
    for irad = 1:nrad
        rr = (irad-0.5)*dh + Ra;
        ll = pi*2*rr;
        ncir = floor(ll / dh);
        tt = (0:ncir-1) .* (pi*2/ncir);
        
        partnum = partnum + ncir;
        xpart(end+1:end+ncir,1) = rr .* cos(tt);
        ypart(end+1:end+ncir,1) = rr .* sin(tt);
        volpart(end+1:end+ncir,1) = pi*((rr+dh/2)^2-(rr-dh/2)^2) / ncir;
        % surfflag(end+1:end+ncir,1) = logical(irad==nrad);
    end
    % surfflag = logical(surfflag);
    % surflist = find(surfflag)';
    % nsurf = numel(surflist);
else
    % grid align
    
    for nx = -nrad:nrad
    for ny = -nrad:nrad
        xx = dh * nx;
        yy = dh * ny;
        rr = sqrt(xx^2 + yy^2);
        
        if rr>Ra && rr<Rb
            partnum = partnum + 1;
            xpart(end+1,1) = xx;
            ypart(end+1,1) = yy;
            volpart(end+1,1) = dh^2;
        end
    end
    end    
end


for ipart = 1:partnum
    isurf = ipart;
% for isurf = 1:nsurf
    % ipart = surflist(isurf);
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


surfflag = (npart<0.95*nzero);
surfflag(xpart.^2+ypart.^2<(Ra+Rb)/2) = 0;
surflist = find(surfflag)';
nsurf = numel(surflist);


if 1
    figure;
    plot(xpart,ypart,'o', xpart(surfflag),ypart(surfflag),'x');
    axis equal;
    hold on;
    rectangle('Position',[-Rb,-Rb,Rb*2,Rb*2],'Curvature',[1 1]);
    hold off;
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



qxfunc = @(x,y) -x ./ sqrt(x.^2+y.^2);
qyfunc = @(x,y) -y ./ sqrt(x.^2+y.^2);
qdivfunc = @(x,y) -1 ./ sqrt(x.^2+y.^2);

spart = [];
if 1
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
        
        % divinner = ndim/nzero * sum(qqs(~surfflag(ineigh)));
        divinner = 2*ndim/nzero * sum(qqs(~surfflag(ineigh)));
        % divinner = ndim/npart(ipart) * sum(qqs(~surfflag(ineigh)));
        % divinner = 2*ndim/nzero * sum(qqs(:));
        
        divqext = qdivfunc(xpart(ipart),ypart(ipart));
        
        spart(isurf,1) = -divqext + divinner;
        % spart(isurf,1) = -divqext;
    end
end
svolpart = spart .* volpart(surfflag);
sum(svolpart)

spart = [];
if 1
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
        
        % spart(isurf,1) = -ndim/nzero * sum(qqs);
        spart(isurf,1) = -2*ndim/nzero * sum(qqs);
        % spart(isurf,1) = -ndim/npart(ipart) * sum(qqs);
    end
end
svolpart = spart .* volpart(surfflag);
sum(svolpart)

if 0
    spart = [];
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
    svolpart = spart .* volpart(surfflag);
    sum(svolpart)
end


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
    if 1
        
        qxs = qxfunc(xpart(ineigh),ypart(ineigh));
        qys = qyfunc(xpart(ineigh),ypart(ineigh));
        
        divq = ((qxs-qxi).*ineighrx + (qys-qyi).*ineighry) ./ ineighr2 .* ineighw;
        % divq = ndim/nzero * sum(divq);
        divq = ndim/npart(ipart) * sum(divq);
    else
        % ineighx = (xpart(ineigh) + xpart(ipart)) / 2;
        % ineighy = (ypart(ineigh) + ypart(ipart)) / 2;
        % qxs = qxfunc(ineighx,ineighy);
        % qys = qyfunc(ineighx,ineighy);
        ineighx = xpart(ineigh);
        ineighy = ypart(ineigh);
        qxs = (qxfunc(ineighx,ineighy) + qxi)/2;
        qys = (qyfunc(ineighx,ineighy) + qyi)/2;
        
        
        divq = ((qxs).*ineighrx + (qys).*ineighry) ./ ineighr2 .* ineighw;
        % divq = 2*ndim/nzero * sum(divq);
        divq = 2*ndim/npart(ipart) * sum(divq);
    end
    [divq, qdivfunc(xpart(ipart),ypart(ipart))]
end













