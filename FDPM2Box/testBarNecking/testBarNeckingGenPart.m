
if prob_type==2 || prob_type==1
    % axisymmetric
    % simply generate a rectangular block
    
    nrad1 = nrad + 1;
    % nlen1 = nlen + 1;
    nlen1 = nmid + 1;
    
    loaded_nodes = [];
    loaded_force = [];
    fixed_nodes_y = [];
    fixed_disp_y = [];
    fixed_nodes_x = [];
    fixed_disp_x = [];
    
    numEdges = 0;
    edgePair = [];
    
    for jlen = 1:nlen1
        % zz = -Lmid + (jlen-1)*(Lmid/nmid); % z coord
        zz = 0 + (jlen-1)*(Lmid/nmid); % z coord
        
        if zz >= 0
            rz = R0*(zz/Lmid) + Rmid*(1-zz/Lmid);
        else
            rz = R0*(-zz/Lmid) + Rmid*(1+zz/Lmid);
        end
        
        for irad = 1:nrad1
            rr = rz/nrad * (irad-1);
            
            numNodes = numNodes + 1;
            
            nodeX(numNodes,1) = rr;
            nodeY(numNodes,1) = zz;
            
            if jlen==1 % bottom end
                fixed_nodes_y(end+1,1) = numNodes;
                fixed_disp_y(end+1,1) = 0;
                % fixed_nodes_x(end+1,1) = numNodes;
                % fixed_disp_x(end+1,1) = 0;
            end
            if jlen==nlen1 % top end
                fixed_nodes_y(end+1,1) = numNodes;
                fixed_disp_y(end+1,1) = Ldisp;
                fixed_nodes_x(end+1,1) = numNodes;
                fixed_disp_x(end+1,1) = 0;
            end
            if irad == 1 % on axis
                % if jlen~=1 && jlen~=nlen1
                if jlen~=nlen1
                % if 1
                    fixed_nodes_x(end+1,1) = numNodes;
                    fixed_disp_x(end+1,1) = 0;
                end
            end
            
            % 
            if jlen==1 && irad<nrad1
                numEdges = numEdges + 1;
                edgePair(numEdges,:) = [numNodes, numNodes+1];
            end
            if irad==nrad1 && jlen<nlen1
                numEdges = numEdges + 1;
                edgePair(numEdges,:) = [numNodes, numNodes+nrad1];
            end
            if jlen==nlen1 && irad<nrad1
                numEdges = numEdges + 1;
                edgePair(numEdges,:) = [numNodes+1, numNodes];
            end
            if irad==1 && jlen<nlen1
                numEdges = numEdges + 1;
                edgePair(numEdges,:) = [numNodes+nrad1, numNodes];
            end
        end
    end
    
    % setup Dirichlet BC
    tmp = fdpmNodeDof(fixed_nodes_y,'y');
    dispBCDofs = [dispBCDofs; tmp];
    dispBCVals = [dispBCVals; fixed_disp_y];
    tmp = fdpmNodeDof(fixed_nodes_x,'x');
    dispBCDofs = [dispBCDofs; tmp];
    dispBCVals = [dispBCVals; fixed_disp_x];
    
    % setup Loading
    tracBCDofs = fdpmNodeDof(loaded_nodes,'xy');
    tracBCVals = reshape(loaded_force, [],1);
end



