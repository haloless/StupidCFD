% Build particle connectivity
% 
%

clear conn;

for i = 1:numNodes
	if 1
        [neigh,rs,rx,ry,cutoff] = fdpmNeighborhood(i,nodePos,re);
        neigh2 = [neigh; i]; % particle included in the neighbourhood
        
        % weight function
        w = fdpmWeightFunc(rs(neigh),cutoff);
        
        % derivative shape function
        dN = fdpmShapeDer(rx(neigh),ry(neigh),w);
        
        % finite-increment-gradient correction to 1st-derivatives
        if exist('fig_stab_alpha') && (fig_stab_alpha > 0)
            % hX = fig_stab_alpha*nodeSize(i,1);
            % hY = fig_stab_alpha*nodeSize(i,2);
            hX = fig_stab_alpha*h0;
            hY = fig_stab_alpha*h0;
            
            [dNX,dNY] = fdpmShapeDerFIGStab(dN(1,:),dN(2,:),dN(3,:),dN(4,:),dN(5,:), hX,hY);
            dN(1,:) = dNX;
            dN(2,:) = dNY;
            % [dN(1,:),dN(2,:)] = fdpmShapeDerFIGStab(dN(1,:),dN(2,:),dN(3,:),dN(4,:),dN(5,:), hX,hY);
        end
    end
	    
	% save connection
    conn(i).numNeigh = numel(neigh2);
	conn(i).neigh2 = neigh2'; % neighbor list is a row vector
    conn(i).W = [w(:)', 0]; % weights, the last self weight is zero
	conn(i).dNX = dN(1,:);
	conn(i).dNY = dN(2,:);
	conn(i).dNXX = dN(3,:);
	conn(i).dNXY = dN(4,:);
	conn(i).dNYY = dN(5,:);
end

if (exist('fig_stab_alpha') && fig_stab_alpha>0 && false)
    
    hvec = fig_stab_alpha * [h0;h0];
    
    for i = 1:numNodes
        
        nn = conn(i).numNeigh;
        neigh2 = conn(i).neigh2;
        dNX = conn(i).dNX;
        dNY = conn(i).dNY;
        
        for ii = 1:conn(i).numNeigh
            j = conn(i).neigh2(ii);
            Hij = [conn(i).dNXX(ii),conn(i).dNXY(ii);conn(i).dNXY(ii),conn(i).dNYY(ii)];
            hij = Hij * hvec;
            dNX(ii) = dNX(ii) + hij(1);
            dNY(ii) = dNY(ii) + hij(2);
        end
        
        for ii = 1:conn(i).numNeigh
            
            k = conn(i).neigh2(ii);
            
            gik = [ conn(i).dNX(ii); conn(i).dNY(ii) ];
            
            for kk = 1:conn(k).numNeigh
                
                j = conn(k).neigh2(kk);
                
                gkj = [ conn(k).dNX(kk); conn(k).dNY(kk) ];
                
                hij = gkj * gik' * hvec;
                
                indj = find(neigh2==j, 1);
                if indj
                    dNX(indj) = dNX(indj) - hij(1);
                    dNY(indj) = dNY(indj) - hij(2);
                else
                    nn = nn + 1;
                    neigh2(end+1) = j;
                    dNX(end+1) = -hij(1);
                    dNY(end+1) = -hij(2);
                end
            end
        end
        
        conn2(i).numNeigh = nn;
        conn2(i).neigh2 = neigh2;
        conn2(i).dNX = dNX;
        conn2(i).dNY = dNY;
    end
    
    clear conn;
    conn = conn2;
    clear conn2;
end










