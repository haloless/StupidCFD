
clear all;

xlo = 0.0;
xhi = 1.0;

nlevel = 18;
np_per_leaf = 4;
np = 2^(nlevel-1) * np_per_leaf;
disp(['np=',int2str(np)]);

%%
%% generate random points
%%
if (1)
    % reset seed
    rng('default');
end
xsrc = xlo + (xhi-xlo).*rand(np,1);
qsrc = rand(np,1);

%%
%% create tree
%%
leveldh = zeros(nlevel,1);
levelnum = zeros(nlevel,1);
levelbeg = zeros(nlevel,1);
levelend = zeros(nlevel,1);
for ilevel = 1:nlevel
    num = LevelNodeNum1D(ilevel);
    range = LevelNodeRange1D(ilevel);
    levelnum(ilevel) = num;
    leveldh(ilevel) = (xhi-xlo) / num;
    levelbeg(ilevel) = range(1);
    levelend(ilevel) = range(end);
end

nnode = 2^nlevel - 1;
xbeg = zeros(nnode,1);
xend = zeros(nnode,1);
xcen = zeros(nnode,1);
nodeside = zeros(nnode,1);  % left=0, right=1
% root
xbeg(1) = xlo;
xend(1) = xhi;
xcen(1) = 0.5*(xlo+xhi);
nodeside(1) = -1;
%
for ilevel = 1:nlevel-1
    for iparent = LevelNodeRange1D(ilevel)
        ichild1 = iparent*2;
        xbeg(ichild1) = xbeg(iparent);
        xend(ichild1) = xcen(iparent);
        xcen(ichild1) = 0.5 * (xbeg(ichild1)+xend(ichild1));
        nodeside(ichild1) = 0;
        ichild2 = iparent*2+1;
        xbeg(ichild2) = xcen(iparent);
        xend(ichild2) = xend(iparent);
        xcen(ichild2) = 0.5 * (xbeg(ichild2)+xend(ichild2));
        nodeside(ichild2) = 1;
    end
end

% register points
regsrc = zeros(np,nlevel);
nodesrc = cell(nnode,1);
for ilevel = 1:nlevel
    num = levelnum(ilevel);
    dh = leveldh(ilevel);
    for ip = 1:np
        inode = floor((xsrc(ip)-xlo)/dh);
        inode = min(inode,num-1) + 2^(ilevel-1);
        regsrc(ip,ilevel) = inode;
        
        % node-particle mapping is hold only on leaf level
        if ilevel == nlevel
            nodesrc{inode}(end+1) = ip;
        end
    end
end



epsilon = 1.0e-6;
mcut = ceil(-log2(epsilon));
disp(['eps=',num2str(epsilon),', mcut=',int2str(mcut)]);


wgt = zeros(nnode,mcut);
disp(['Compute weight']);
tic;
for imom = 1:mcut
    m = imom - 1;
    for ilevel = 1:nlevel
        for ip = 1:np
            inode = regsrc(ip,ilevel);
            am = (-1)^m * (m+1) * (xsrc(ip)-xcen(inode))^m;
            % TODO multiply q here
            wgt(inode,imom) = wgt(inode,imom) + am*qsrc(ip);
        end
    end
end
toc;

% ydst = rand(1,1);
ydst = 0.73;
disp(['ydst=',num2str(ydst)]);

ystclfar = [];
ystclnear = [];
ysumfar = 0;
ysumnear = 0;
disp(['Compute FMM']);
tic;
for ilevel = 1:nlevel
    
    iloc = floor((ydst-xlo)/leveldh(ilevel));
    iloc = max(0,min(iloc,levelnum(ilevel)-1)) + 2^(ilevel-1);
    disp(['level=',int2str(ilevel), '; loc=',int2str(iloc)]);
    
    if ilevel<nlevel && ilevel>1
        % on coarse level
        % find interaction list
        ipar = floor(iloc/2);
        nearlo = max(ipar-1,levelbeg(ilevel-1));
        nearhi = min(ipar+1,levelend(ilevel-1));
        actlo = nearlo*2;
        acthi = nearhi*2+1;
        active = actlo:acthi;
        interact = active(active<iloc-1 | active>iloc+1);
        % disp(['far interact=',int2str(interact)]);
        ystclfar = [ystclfar, interact];
        
        % add far-field
%         for imom = 1:mcut
%             m = imom - 1;
%             for jnode = interact
%                 wm = wgt(jnode,imom);
%                 rr = xcen(jnode) - ydst;
%                 sm = 1.0/rr^2 / rr^m;
%                 ysumfar = ysumfar + wm*sm;
%             end
%         end
    elseif ilevel==nlevel
        % on finest level
        % add near-field
        
        if nodeside(iloc) == 0
            nearlo = iloc - 2;
            nearhi = iloc + 3;
        else
            nearlo = iloc - 3;
            nearhi = iloc + 2;
        end
        nearlo = max(nearlo,levelbeg(ilevel));
        nearhi = min(nearhi,levelend(ilevel));
        active = nearlo:nearhi;
        
        interact = active(active<iloc-1 | active>iloc+1);
        % disp(['far interact=',int2str(interact)]);
        ystclfar = [ystclfar, interact];
        
        interact = active(active>=iloc-1 & active<=iloc+1);
        % disp(['near interact=',int2str(interact)]);
        ystclnear = [ystclnear, interact];
        
%         for jnode = interact
%             for jp = find(regsrc(:,ilevel)==jnode)' 
%                 %[jnode,jp]
%                 ysumnear = ysumnear + 1.0/(xsrc(jp)-ydst)^2 * qsrc(jp);
%             end
%         end
    end
    
end
for jnode = ystclfar
    for imom = 1:mcut
        m = imom - 1;
        wm = wgt(jnode,imom);
        rr = xcen(jnode) - ydst;
        sm = 1.0/rr^2 / rr^m;
        ysumfar = ysumfar + wm*sm;
    end
end
for jnode = ystclnear
    for jp = nodesrc{jnode} % find(regsrc(:,ilevel)==jnode)' 
        %[jnode,jp]
        ysumnear = ysumnear + 1.0/(xsrc(jp)-ydst)^2 * qsrc(jp);
    end
end
ysum = ysumfar + ysumnear;
toc;
disp(['FMM: ysum=',num2str(ysum), ', far=',num2str(ysumfar), ', near=',num2str(ysumnear)]);


% direct sum check
disp(['Compute Direct']);
ysum0 = 0;
tic;
for jp = 1:np
    ysum0 = ysum0 + qsrc(jp)/(ydst-xsrc(jp))^2;
end
toc;
yerr = abs(ysum-ysum0) / ysum0;
disp(['Direct: ysum=',num2str(ysum0), ', far=',num2str(ysum0-ysumnear), ', near=',num2str(ysumnear)]);
disp(['err=',num2str(yerr)]);


if (0)
    figure;
    hold on;
    for ilevel = 1:nlevel
        num = LevelNodeNum1D(ilevel);
        range = LevelNodeRange1D(ilevel);
        
        % plot grid
        xx = xcen(range);
        ll = ilevel.*ones(size(xx));
        plot(xx,ll,'o');
        xx = reshape([xbeg(range),xend(range)]',num*2,1);
        ll = ilevel.*ones(size(xx));
        plot(xx,ll,'+-');
        for inode = range
            text(xcen(inode),ilevel,int2str(inode), ...
            'FontSize',14,'Color','r','VerticalAlignment','bottom','HorizontalAlignment','center');
        end
        
        % plot source
        if (0)
        pp = cell(3,num);
        for inode = 1:num
            pp{1,inode} = xsrc(find(regsrc(:,ilevel)==range(inode)));
            pp{2,inode} = ilevel.*ones(size(pp{1,inode}));
            pp{3,inode} = '.';
        end
        plot(pp{:});
        end
        
        % plot target
        plot(ydst,ilevel,'xr');
    end
    
    hold off;
    axis([xlo-0.1,xhi+0.1,0,nlevel+1]);
end

if (0)
    figure;
    hold on;
    % plot active stencil
    nstcl = size(ystclfar(:),1);
    for istcl = ystclfar
        ilevel = floor(log2(istcl))+1;
        xx = [xbeg(istcl),xend(istcl)];
        ll = ilevel.*ones(size(xx)); 
        ll(:) = 1;
        plot(xx,ll,'.-b');
    end
    nstcl = size(ystclnear(:),1);
    for istcl = ystclnear
        ilevel = floor(log2(istcl))+1;
        xx = [xbeg(istcl),xend(istcl)];
        ll = ilevel.*ones(size(xx));
        ll(:) = 1;
        plot(xx,ll,'.-r');
    end
    plot(ydst,nlevel,'xk');
    hold off;
    axis([xlo-0.1,xhi+0.1,0,nlevel+1]);
end











