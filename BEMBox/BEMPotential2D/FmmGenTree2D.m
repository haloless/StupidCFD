
ielem = zeros(n,1);
itree = zeros(ncellmx,1);
loct = zeros(ncellmx,1);
numt = zeros(ncellmx,1);
ifath = zeros(ncellmx,1);
leafflag = zeros(ncellmx,1);

% NOTE level have zero index
level = zeros(levmx+1,1); 
levelndiv = zeros(levmx+1,1);
leveldx = zeros(levmx+1,1);
leveldy = zeros(levmx+1,1);

% initialize to map all elements
ielem(:) = 1:n;

% the root level cell
itree(1) = 0;
level(1) = 1;
level(2) = 2;
loct(1) = 1;
ifath(1) = 1;
numt(1) = n;
leafflag(1) = 0;

levelndiv(1) = 1;
leveldx(1) = (xmax-xmin);
leveldy(1) = (ymax-ymin);

ndivx = 1;
lowlev = 1;
nleaf = 0;

for lev = 1:levmx
    levp = lev - 1;
    levn = lev + 1;
    level(levn+1) = level(lev+1);
    if level(lev+1) == level(levp+1); break; end
    
    % refine cell size
    ndivxp = ndivx;
    ndivx = 2*ndivxp;
    dxp = (xmax-xmin) / ndivxp;
    dyp = (ymax-ymin) / ndivxp;
    
    levelndiv(lev+1) = ndivx;
    leveldx(lev+1) = dxp/2;
    leveldy(lev+1) = dyp/2;
    
    for inp = level(levp+1):level(lev+1)-1
        itrp = itree(inp);
        
        % the parent cell can be further refined
        if numt(inp)>maxl || (lev<=2 && numt(inp)>0)
            % cell index in the parent level
            itrpx = mod(itrp,ndivxp);
            itrpy = floor(itrp/ndivxp);
            
            % center point of parent cell
            % use this to divide parent
            xsep = xmin + (itrpx+0.5)*dxp;
            ysep = ymin + (itrpy+0.5)*dyp;
            
            % do bisect to separate the parent cell
            % separate by y
            [ielem,nsepy] = FmmBisect(x,ielem,loct(inp), numt(inp), ysep, 2);
            % separate lower part by x
            [ielem,nsepx1] = FmmBisect(x,ielem,loct(inp), nsepy-1, xsep, 1);
            % separate upper part by x
            [ielem,nsepx2] = FmmBisect(x,ielem,loct(inp)+nsepy-1, numt(inp)-(nsepy-1), xsep,1);
            
            % make four child cells
            nwk = [ nsepx1-1, nsepy-nsepx1, nsepx2-1, numt(inp)-nsepy-nsepx2+2 ];
            locc = loct(inp);
            
            for icldy = 0:1
            for icldx = 0:1
                % local child index 1~4
                icld = icldx + icldy*2 + 1;
                if nwk(icld) > 0
                    nrel = level(levn+1);
                    if nrel > ncellmx
                        error('ncell overflow');
                    end
                    
                    % cell index in the child level
                    itrx = itrpx*2 + icldx;
                    itry = itrpy*2 + icldy;
                    itree(nrel) = itry*ndivx + itrx;
                    
                    loct(nrel) = locc;
                    numt(nrel) = nwk(icld);
                    ifath(nrel) = inp;
                    lowlev = lev;
                    
                    % child cell is leaf
                    if lev>1 && (numt(nrel)<=maxl || lev==levmx)
                        nleaf = nleaf + 1;
                        if nleaf > nleafmx
                            error('nleaf overflow');
                        end
                        leafflag(nrel) = 1;
                    else
                        leafflag(nrel) = 0;
                    end
                    
                    level(levn+1) = nrel + 1;
                    locc = locc + nwk(icld);
                end
            end
            end
        end
    end
end

if (1)
    disp(['number of tree levels = ', int2str(lowlev)]);
    disp(['number of tree leaves = ', int2str(nleaf)]);
    disp(['number of tree cells = ', int2str(nrel)]);
end
if (0)
    disp('quadtree structure');
    fprintf(1, '%8s %8s %8s %8s %8s %8s\n', 'cellno','itree','loct','numt','ifath','ielem');
    I = 1:nrel;
    fprintf(1, '%8d %8d %8d %8d %8d %8d\n', [I',itree(I),loct(I),numt(I),ifath(I),ielem(I)]');
end
if (1)
    % plot all levels
    figure;
    hold on;
    
    for lev = 0:levmx
        thislev = lev+1;
        nextlev = lev+2;
        if level(thislev)==level(nextlev); break; end
        
        ndiv = 2^lev;
        dx = (xmax-xmin) / ndiv;
        dy = (ymax-ymin) / ndiv;
        
        for in = level(thislev):level(nextlev)-1
            itr = itree(in);
            itrx = mod(itr,ndiv);
            itry = floor(itr/ndiv);
            x0 = xmin + itrx*dx;
            y0 = ymin + itry*dy;
            x1 = x0 + dx;
            y1 = y0 + dy;
            plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],'.-');
            text((x0+x1)/2,(y0+y1)/2,int2str(in));
        end
    end
    
    plot(y(1,:),y(2,:),'.-');
    % plot(x(1,:),x(2,:),'o');
    
    hold off;
    axis equal;
    axis([xmin-xlen*0.05, xmax+xlen*0.05 ymin-ylen*0.05 ymax+ylen*0.05]);
end
if (0)
    % plot each level
    for lev = 0:levmx
        thislev = lev+1;
        nextlev = lev+2;
        if level(thislev)==level(nextlev); break; end
        
        figure;
        hold on;
        
        ndiv = 2^lev;
        dx = (xmax-xmin) / ndiv;
        dy = (ymax-ymin) / ndiv;
        
        for in = level(thislev):level(nextlev)-1
            itr = itree(in);
            itrx = mod(itr,ndiv);
            itry = floor(itr/ndiv);
            x0 = xmin + itrx*dx;
            y0 = ymin + itry*dy;
            x1 = x0 + dx;
            y1 = y0 + dy;
            % draw cell
            plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],'.-');
            text((x0+x1)/2,(y0+y1)/2,int2str(in));
            
            iloc = loct(in):loct(in)+numt(in)-1;
            % plot(x(1,elem(iloc)),x(2,elem(iloc)),'x');
            for il = iloc
                ip = ielem(il);
                text(x(1,ip),x(2,ip), int2str(in));
            end
        end
        
        hold off;
        axis equal;
        axis([xmin-xlen*0.05, xmax+xlen*0.05 ymin-ylen*0.05 ymax+ylen*0.05]);
        title(['level ',int2str(lev)]);
    end
end








