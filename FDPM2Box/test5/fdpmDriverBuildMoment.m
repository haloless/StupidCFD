%%fdpmDriverBuildMoment
% Build moment

nodeMomX = zeros(numNodes,1);
nodeMomY = zeros(numNodes,1);
nodeMomXX = zeros(numNodes,1);
nodeMomYY = zeros(numNodes,1);
nodeMomXY = zeros(numNodes,1);


%
% gaussn = 3;
gaussn = 5;
% gaussn = 7;
[gaussx,gaussw] = gauss1dPoint(gaussn);

nodeFlag = zeros(numNodes, 1);
nodeFlag(trep.freeBoundary()) = 1;


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
        
        for k = 1:gaussn
            xk = gp(k,1);
            yk = gp(k,2);
            wk = len * gaussw(k) / sum(gaussw);
            
            wk1 = wk;
            if use_rz_coord, wk1 = wk * nodeX(i1); end
            nodeMomX(i1) = nodeMomX(i1) + 1/2*(xk-nodeX(i1))^2 * nvec(1) * wk1;
            nodeMomY(i1) = nodeMomY(i1) + 1/2*(yk-nodeY(i1))^2 * nvec(2) * wk1;
            nodeMomXX(i1) = nodeMomXX(i1) + 1/3*(xk-nodeX(i1))^3 * nvec(1) * wk1;
            nodeMomYY(i1) = nodeMomYY(i1) + 1/3*(yk-nodeY(i1))^3 * nvec(2) * wk1;
            nodeMomXY(i1) = nodeMomXY(i1) + 1/2*(xk-nodeX(i1))^2*(yk-nodeY(i1)) * nvec(1)*wk1;
            
            wk2 = wk;
            if use_rz_coord, wk2 = wk * nodeX(i2); end
            nodeMomX(i2) = nodeMomX(i2) - 1/2*(xk-nodeX(i2))^2 * nvec(1) * wk2;
            nodeMomY(i2) = nodeMomY(i2) - 1/2*(yk-nodeY(i2))^2 * nvec(2) * wk2;
            nodeMomXX(i2) = nodeMomXX(i2) - 1/3*(xk-nodeX(i2))^3 * nvec(1) * wk2;
            nodeMomYY(i2) = nodeMomYY(i2) - 1/3*(yk-nodeY(i2))^3 * nvec(2) * wk2;
            nodeMomXY(i2) = nodeMomXY(i2) - 1/2*(xk-nodeX(i2))^2*(yk-nodeY(i2)) * nvec(1)*wk2;
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
                    
                    if use_rz_coord, wk = wk * nodeX(ind); end
                    
                    nodeMomX(ind) = nodeMomX(ind) + 1/2*(xk-nodeX(ind))^2 * nvec(1) * wk;
                    nodeMomY(ind) = nodeMomY(ind) + 1/2*(yk-nodeY(ind))^2 * nvec(2) * wk;
                    nodeMomXX(ind) = nodeMomXX(ind) + 1/3*(xk-nodeX(ind))^3 * nvec(1) * wk;
                    nodeMomYY(ind) = nodeMomYY(ind) + 1/3*(yk-nodeY(ind))^3 * nvec(2) * wk;
                    nodeMomXY(ind) = nodeMomXY(ind) + 1/2*(xk-nodeX(ind))^2*(yk-nodeY(ind)) * nvec(1)*wk;
                end
            end
        end
    end
    
end



