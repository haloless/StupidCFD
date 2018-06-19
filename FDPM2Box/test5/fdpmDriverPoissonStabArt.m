
cart = -0.05;

for i = 1:numNodes
    % nneigh = conn(i).numNeigh;
    ineigh = conn(i).neigh2;
    ivol = nodeVol(i);
    
    
    %
    for j = ineigh
        jneigh = conn(j).neigh2;
        jvol = nodeVol(j);
        
        rij = nodePos(j,:) - nodePos(i,:);
        wij = 4 - norm(rij)^2 / re^2;
        
        cij = cart;
        
        %
        Kassem.SpMatAssemBlockWithDof([1,-1].*wij.*cij, i,[j,i]);
        
        %
        Kassem.SpMatAssemBlockWithDof(-rij * [conn(i).dNX; conn(i).dNY] .* wij.*cij, i,ineigh);
        
        %
        Kassem.SpMatAssemBlockWithDof(-rij * [conn(j).dNX; conn(j).dNY] .* wij.*cij, i,jneigh);
    end
end





