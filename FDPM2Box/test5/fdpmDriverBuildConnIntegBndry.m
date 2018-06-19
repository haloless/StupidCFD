%Build connection

clear conn;

%
for i = 1:numNodes
    conn(i).numNeigh = 0;
    conn(i).avol = nodeVol(i);
    
    conn(i).tmp = zeros(2,numNodes);
    conn(i).neigh = zeros(numNodes,1);
    
    % conn(i).neigh1 = find(sqrt((nodeX(i)-nodeX).^2 + (nodeY(i)-nodeY).^2) < re*2);
end

%
gaussn = 7;
[gaussx,gaussw] = gauss1dPoint(gaussn);

nodeFlag = zeros(numNodes, 1);
nodeFlag(deltri.freeBoundary()) = 1;

figure;
triplot(tri, nodeX, nodeY);
axis equal;
hold on;
plot(nodeX,nodeY,'.k');

%
for i = 1:numTri
    
    % A has the max angle
    [ia,ib,ic, la,lb,lc, ta,tb,tc] = triRotateMaxAngle(tri, nodePos, i);
    
    pa = nodePos(ia,:);
    pb = nodePos(ib,:);
    pc = nodePos(ic,:);
    
    % internal boundary
    if ta >= pi/2 - 1.0e-5
        pout = 0.5 * (pb + pc);
        pairs = [ ia,ib; ia,ic ];
    else
        %
        pout = triCircumCenter(pa(1),pa(2), pb(1),pb(2), pc(1),pc(2));
        pairs = [ ia,ib; ib,ic; ic,ia ];
    end
    npairs = size(pairs,1);
    for ii = 1:npairs
        i1 = pairs(ii,1);
        i2 = pairs(ii,2);
        
        pmid = 0.5 * (nodePos(i1,:) + nodePos(i2,:));
        len = norm(pout - pmid);
        
        nvec = nodePos(i2,:) - nodePos(i1,:);
        nvec = nvec - dot(nvec,pout-pmid)./len.*(pout-pmid)./len;
        nvec = nvec ./ norm(nvec);
        
        gp = ((gaussx+1)./2) * pout + ((1-gaussx)./2) * pmid;
        
        if npairs == 2
            plot(gp(:,1), gp(:,2), '.g');
        else
            plot(gp(:,1), gp(:,2), '.r');
        end
        
        for k = 1:gaussn
            xk = gp(k,1);
            yk = gp(k,2);
            wk = len * gaussw(k) / sum(gaussw);
            
            ok = find(sqrt((xk-nodeX).^2+(yk-nodeY).^2) < re);
            
            [N,Nx,Ny] = mlsRegularShape2D([xk,yk], nodePos(ok,1),nodePos(ok,2), re);
            
            
            conn(i1).neigh(ok) = 1;
            conn(i1).tmp(:,ok) = conn(i1).tmp(:,ok) + nvec' * N .* wk;
            
            conn(i2).neigh(ok) = 1;
            conn(i2).tmp(:,ok) = conn(i2).tmp(:,ok) - nvec' * N .* wk;
        end
    end
    
    
    % free boundary
    for ii = 1:3
        if ii == 1
            i1 = ia; i2 = ib; i3 = ic;
        elseif ii == 2
            i1 = ib; i2 = ic; i3 = ia;
        else
            i1 = ic; i2 = ia; i3 = ib;
        end
        if nodeFlag(i1) && nodeFlag(i2)
            pmid = 0.5 * (nodePos(i1,:) + nodePos(i2,:));
            vdif = nodePos(i2,:) - nodePos(i1,:);
            vlen = norm(vdif);
            
            nvec = vdif ./ vlen;
            nvec = [ nvec(2), -nvec(1) ];
            
            for jj = 1:2
                if jj == 1
                    ind = i1;
                    pbeg = nodePos(i1,:);
                    pend = pmid;
                else
                    ind = i2;
                    pbeg = pmid;
                    pend = nodePos(i2,:);
                end
                
                gp = ((gaussx+1)./2) * pend + ((1-gaussx)./2) * pbeg;
                len = vlen / 2;
                
                for k = 1:gaussn
                    xk = gp(k,1);
                    yk = gp(k,2);
                    wk = len * gaussw(k) / sum(gaussw);
                    
                    ok = find(sqrt((xk-nodeX).^2+(yk-nodeY).^2) < re);
                    
                    [N,Nx,Ny] = mlsRegularShape2D([xk,yk], nodePos(ok,1),nodePos(ok,2), re);
                    
                    conn(ind).neigh(ok) = 1;
                    conn(ind).tmp(:,ok) = conn(ind).tmp(:,ok) + nvec' * N .* wk;
                    
                end
            end
        end
    end
    
end


% clean up
for i = 1:numNodes
    conn(i).neigh2 = find(conn(i).neigh).';
    if 1 % swap self to the last position
        ii = find(conn(i).neigh2 == i);
        if ii
            % tmp = conn(i).neigh2(end); conn(i).neigh2(end) = conn(i).neigh2(ii); conn(i).neigh2(ii) = tmp;
            % tmp = conn(i).dNX(end); conn(i).dNX(end) = conn(i).dNX(ii); conn(i).dNX(ii) = tmp;
            % tmp = conn(i).dNY(end); conn(i).dNY(end) = conn(i).dNY(ii); conn(i).dNY(ii) = tmp;
            % tmp = conn(i).N(end); conn(i).N(end) = conn(i).N(ii); conn(i).N(ii) = tmp;
            conn(i).neigh2(ii) = [];
            conn(i).neigh2(end+1) = i;
        end
    end
    
    conn(i).numNeigh = length(conn(i).neigh2);
    conn(i).dNX = conn(i).tmp(1, conn(i).neigh2) ./ conn(i).avol;
    conn(i).dNY = conn(i).tmp(2, conn(i).neigh2) ./ conn(i).avol;
    
    conn(i).neigh = [];
    conn(i).tmp = [];
    
    % nodeVol(i) = conn(i).asum;
end




