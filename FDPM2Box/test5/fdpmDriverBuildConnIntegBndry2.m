%fdpmDriverBuildConnIntegBndry2: Build connection by integration on divided triangle boundaries.
%
% Output
% - nodeVol


%% NOTE on axisymmetric RZ coord
%


clear conn;

if ~exist('use_auto_re','var')
    use_auto_re = 1; % automatically adjust re by triangle size
end

show_tri_div = 1; % plot sub-divided triangles

use_rz_coord = 0; % current coordinate system
if exist('prob_type','var') && prob_type==2 
    % in RZ coord, the integral volume will be corrected
    use_rz_coord = 1;
end

% set volume
nodeVol = nodeArea;

%
for i = 1:numNodes
    conn(i).numNeigh = 0;
    
    % the integral volume
    conn(i).avol = nodeVol(i); 
    if use_rz_coord
        conn(i).avol2 = 0;
    end
    
    % temporary buffers
    conn(i).tmp = zeros(2,numNodes);
    conn(i).tmp2 = zeros(3,numNodes);
    
    conn(i).neigh = zeros(numNodes,1);
end

%
% gaussn = 3;
gaussn = 5;
% gaussn = 7;
% gaussn = 11;
[gaussx,gaussw] = gauss1dPoint(gaussn);

% tag nodes on free boundaries
nodeFlag = zeros(numNodes, 1);
nodeFlag(trep.freeBoundary()) = 1;

if show_tri_div % draw divided triangles
    figure;
    triplot(tri, nodeX, nodeY);
    axis equal;
    hold on;
    plot(nodeX,nodeY,'.k');
end

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
        
        if show_tri_div
            if npairs == 2
                plot(gp(:,1), gp(:,2), '.g');
            else
                plot(gp(:,1), gp(:,2), '.r');
            end
        end
        
        if use_auto_re % variable re
            re = dilation * norm(nodePos(i2,:)-nodePos(i1,:));
        end
        
        for k = 1:gaussn
            xk = gp(k,1);
            yk = gp(k,2);
            wk = len * gaussw(k) / sum(gaussw);
            
            ok = find(sqrt((xk-nodeX).^2+(yk-nodeY).^2) < re);
            
            [N,Nx,Ny] = mlsRegularShape2D([xk,yk], nodePos(ok,1),nodePos(ok,2), re);
            
            wk1 = wk;
            % if use_rz_coord, wk1 = wk * xk; end
            wk2 = wk;
            % if use_rz_coord, wk2 = wk * xk; end;
            
            % gradient
            conn(i1).neigh(ok) = 1;
            conn(i1).tmp(:,ok) = conn(i1).tmp(:,ok) + nvec' * N .* wk1;
            
            conn(i2).neigh(ok) = 1;
            conn(i2).tmp(:,ok) = conn(i2).tmp(:,ok) - nvec' * N .* wk1;
            
            % hessian
            conn(i1).tmp2(1,ok) = conn(i1).tmp2(1,ok) + nvec(1).*Nx.*wk2;
            conn(i1).tmp2(2,ok) = conn(i1).tmp2(2,ok) + nvec(2).*Ny.*wk2;
            conn(i1).tmp2(3,ok) = conn(i1).tmp2(3,ok) + 0.5 .* (nvec(1).*Ny+nvec(2).*Nx) .* wk2;
            
            conn(i2).tmp2(1,ok) = conn(i2).tmp2(1,ok) - nvec(1).*Nx.*wk2;
            conn(i2).tmp2(2,ok) = conn(i2).tmp2(2,ok) - nvec(2).*Ny.*wk2;
            conn(i2).tmp2(3,ok) = conn(i2).tmp2(3,ok) - 0.5 .* (nvec(1).*Ny+nvec(2).*Nx) .* wk2;
            
            % volume
            if use_rz_coord
                conn(i1).avol2 = conn(i1).avol2 + 0.5*xk^2 * nvec(1) * wk;
                conn(i2).avol2 = conn(i2).avol2 - 0.5*xk^2 * nvec(1) * wk;
            end
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
            
            if use_auto_re
                re = dilation * vlen;
            end
            
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
                    
                    wk1 = wk;
                    % if use_rz_coord, wk1 = wk * xk; end
                    wk2 = wk;
                    % if use_rz_coord, wk2 = wk * xk; end;
                    
                    conn(ind).neigh(ok) = 1;
                    conn(ind).tmp(:,ok) = conn(ind).tmp(:,ok) + nvec' * N .* wk1;
                    
                    conn(ind).tmp2(1,ok) = conn(ind).tmp2(1,ok) + nvec(1).*Nx.*wk2;
                    conn(ind).tmp2(2,ok) = conn(ind).tmp2(2,ok) + nvec(2).*Ny.*wk2;
                    conn(ind).tmp2(3,ok) = conn(ind).tmp2(3,ok) + 0.5 .* (nvec(1).*Ny+nvec(2).*Nx) .* wk2;
                    
                    if use_rz_coord % volume
                        conn(ind).avol2 = conn(ind).avol2 + 0.5*xk^2 * nvec(1) * wk;
                    end
                end
            end
        end
    end
    
end


% clean up
for i = 1:numNodes
    ineigh = find(conn(i).neigh).';
    if 1 % swap self to the last position
        ii = find(ineigh == i);
        if ii
            ineigh(ii) = [];
            ineigh(end+1) = i;
        end
    end
    conn(i).neigh2 = ineigh;
    conn(i).numNeigh = length(ineigh);
    
    if 1 % 1st-order derivative
        conn(i).dNX = conn(i).tmp(1, ineigh) ./ conn(i).avol;
        conn(i).dNY = conn(i).tmp(2, ineigh) ./ conn(i).avol;
    end
    
    if 1 % 2nd-order smoother
        conn(i).dNXX = conn(i).tmp2(1,ineigh) ./ conn(i).avol;
        conn(i).dNYY = conn(i).tmp2(2,ineigh) ./ conn(i).avol;
        conn(i).dNXY = conn(i).tmp2(3,ineigh) ./ conn(i).avol;
    end
    
    % corrected volume
    if use_rz_coord
        nodeVol(i) = conn(i).avol2;
        % nodeVol(i) = nodeVol(i) * nodeX(i);
    end
    
    % remove useless 
    conn(i).neigh = [];
    conn(i).tmp = [];
    conn(i).tmp2 = [];
    
end




