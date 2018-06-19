
% fstab,Kstab
% Estab is the coefficient

if 1
    
    fstab = zeros(numDofs,1);
    Kassem = SpMatAssem(numDofs);
    
    for i = 1:numNodes
        nneigh = conn(i).numNeigh;
        ineigh = conn(i).neigh2;
        wneigh = conn(i).W;
        
        idof = fdpmNodeDof(i,'xy');
        idofx = fdpmNodeDof(i,'x');
        idofy = fdpmNodeDof(i,'y');
        
        ineighdof = fdpmNodeDof(ineigh,'xy');
        ineighdofx = fdpmNodeDof(ineigh,'x');
        ineighdofy = fdpmNodeDof(ineigh,'y');
        
        ix = nodeCurr(1,i);
        iy = nodeCurr(2,i);
        ineighx = nodeCurr(1,ineigh)';
        ineighy = nodeCurr(2,ineigh)';
        
        % DN = [conn(i).dNX; conn(i).dNY];
        
        % Fi = ineighpos * DN';
        
        % dN = Fi' \ DN;
        
        % if prob_type == 1
            % Fi(3,3) = 1.0;
        % elseif prob_type == 2
            % Fi(3,3) = ineighpos(1,end) / nodeCoord(1,i);
        % end
        % detF = det(Fi);
        
        Gi = [conn(i).dNX; conn(i).dNY];
        Hi = [conn(i).dNXX; conn(i).dNXY; conn(i).dNYY];
        
        dX = nodeX(ineigh)-nodeX(i);
        dY = nodeY(ineigh)-nodeY(i);
        Mi = fdpmShapeBasis(dX,dY);
        
        Li = Mi * [Gi;Hi];
        Li(:,end) = Li(:,end) + 1;
        
        %
        Ks1 = wneigh * Li - wneigh;
        Ks1 = Ks1 .* Estab;
        SpMatAssemBlockWithDof(Kassem, Ks1, idofx,ineighdofx);
        SpMatAssemBlockWithDof(Kassem, Ks1, idofy,ineighdofy);
        
        %
        Ks2 = bsxfun(@times, Li, wneigh');
        Ks2 = Ks2 .* (-Estab);
        SpMatAssemBlockWithDof(Kassem, Ks2, ineighdofx,ineighdofx);
        SpMatAssemBlockWithDof(Kassem, Ks2, ineighdofy,ineighdofy);
        
        %
        Ks3 = diag(wneigh);
        Ks3 = Ks3 .* Estab;
        SpMatAssemBlockWithDof(Kassem, Ks3, ineighdofx,ineighdofx);
        SpMatAssemBlockWithDof(Kassem, Ks3, ineighdofy,ineighdofy);
        
        % dx/dX
        gx = Gi * ineighx;
        hx = Hi * ineighx;
        gy = Gi * ineighy;
        hy = Hi * ineighy;
        
        xint = Mi * [gx;hx] + ix;
        yint = Mi * [gy;hy] + iy;
        
        %
        fsx = (xint-ineighx).*wneigh' * Estab;
        fsy = (yint-ineighy).*wneigh' * Estab;
        
        fstab(idofx) = fstab(idofx) + sum(fsx);
        fstab(idofy) = fstab(idofy) + sum(fsy);
        
        %
        fstab(ineighdofx) = fstab(ineighdofx) - fsx;
        fstab(ineighdofy) = fstab(ineighdofy) - fsy;
        
    end


    Kstab = SpMatCreateAndClear(Kassem);

end










