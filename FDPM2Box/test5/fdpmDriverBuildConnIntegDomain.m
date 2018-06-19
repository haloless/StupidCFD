%Build connection by domain integration 
% Closest Gauss points are used.

% prepare connections
clear conn;

for i = 1:numNodes
    conn(i).numNeigh = 0;
    conn(i).asum = 0;
    
    conn(i).tmp = zeros(3,numNodes);
    conn(i).neigh = zeros(numNodes,1);
    
    conn(i).neigh1 = find(sqrt((nodeX(i)-nodeX).^2 + (nodeY(i)-nodeY).^2) < re*2);
end

% gaussord = 7;
gaussord = 11;
[gaussp, gaussw] = gauss2dTriPoint(gaussord);
gaussx = gaussp(:,1);
gaussy = gaussp(:,2);
gaussn = length(gaussw);

ntri = size(tri, 1);
for j = 1:ntri
    i1 = tri(j,1);
    i2 = tri(j,2);
    i3 = tri(j,3);
    x21 = nodeX(i2) - nodeX(i1);
    x31 = nodeX(i3) - nodeX(i1);
    y21 = nodeY(i2) - nodeY(i1);
    y31 = nodeY(i3) - nodeY(i1);
    gpx = nodeX(i1) + x21.*gaussx + x31.*gaussy;
    gpy = nodeY(i1) + y21.*gaussx + y31.*gaussy;
    ga = 0.5 * abs(x21*y31 - y21*x31);
    
    for k = 1:gaussn
        % gauss point
        xk = gpx(k);
        yk = gpy(k);
        pk = [xk,yk];
        wk = ga * gaussw(k) / sum(gaussw);
        if prob_type == 2 % axisymmetric
            wk = wk * 2*pi*xk;
        end
        
        % owner is the closest node
        r1 = norm(pk - nodePos(i1,:));
        r2 = norm(pk - nodePos(i2,:));
        r3 = norm(pk - nodePos(i3,:));
        if r1 > r2
            if r2 > r3
                imin = i3;
            else
                imin = i2;
            end
        else
            if r1 > r3
                imin = i3;
            else
                imin = i1;
            end
        end
        
        % ok = find(sqrt((xk-nodeX).^2+(yk-nodeY).^2) < re);
        ok = find(sqrt((xk-nodeX(conn(imin).neigh1)).^2 + (yk-nodeY(conn(imin).neigh1)).^2) < re);
        ok = conn(imin).neigh1(ok);
        
        [N,Nx,Ny] = mlsRegularShape2D([xk,yk], nodePos(ok,1),nodePos(ok,2), re);
        % vsum2(ok,1) = vsum2(ok,1) + wk .* Nx(:);
        % vsum2(ok,2) = vsum2(ok,2) + wk .* Ny(:);
        
        dN = [Nx; Ny; N];
        % Ke = dN' * dN .* wk;
        
        % SpMatAssemBlockWithDof(Kassem, Ke, ok,ok);
        % fext2(ok) = fext2(ok) + fun_fext(xk,yk)*wk .* N(:);
        
        
        
        conn(imin).neigh(ok) = 1;
        conn(imin).tmp(:,ok) = conn(imin).tmp(:,ok) + dN.*wk;
        conn(imin).asum = conn(imin).asum + wk;
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
    conn(i).dNX = conn(i).tmp(1, conn(i).neigh2) ./ conn(i).asum;
    conn(i).dNY = conn(i).tmp(2, conn(i).neigh2) ./ conn(i).asum;
    conn(i).N = conn(i).tmp(end, conn(i).neigh2) ./ conn(i).asum;
    
    
    conn(i).neigh = [];
    conn(i).tmp = [];
    
    nodeVol(i) = conn(i).asum;
end

