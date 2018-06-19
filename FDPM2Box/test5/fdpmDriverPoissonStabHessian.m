
for i = 1:numNodes
    nneigh = conn(i).numNeigh;
    ineigh = conn(i).neigh2;
    ivol = nodeVol(i);
    
    % finite increment
    eta = 1.0;
    hfic = eta * h0;
    
    % cfic = 0.25;
    cfic = 1.0 / 12;
    
    if 0
        dH = [conn(i).dNXX; conn(i).dNXY; conn(i).dNXY; conn(i).dNYY];
        Kfic = dH' * dH * ivol * cfic*hfic^2;
    else
        Bx = [conn(i).dNXX; conn(i).dNXY];
        By = [conn(i).dNXY; conn(i).dNYY];
        Kfic = nodeMomXX(i)*Bx'*Bx + nodeMomYY(i)*By'*By + nodeMomXY(i)*Bx'*By + nodeMomXY(i)*By'*Bx;
    end
    
    SpMatAssemBlockWithDof(Kassem, Kfic, ineigh,ineigh);
end



