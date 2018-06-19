% Generate particles and determine loadings
%
%

if prob_type == 1
    % plane strain
    % generate a circular slice of the tube
    
    nlen = round(0.5*pi*Ra/h0);
    nlen1 = nlen + 1;
    nrad1 = nrad + 1;
    
    % used to find most inner particles for apply pressure
    loaded_nodes = [];
    loaded_force = [];
    
    % tag inner surface
    numEdges = 0;
    edgePair = [];
    
    for irad = 1:nrad1
        % 1/4 cylinder, constraint on x- and y-axis
        rr = Ra + (irad-1)*h0;
        
        nn = nlen1;
        for ii = 1:nn
            tt = pi/2 / (nn-1) * (ii-1);
            
            numNodes = numNodes + 1;
            nodeX(numNodes,1) = rr * cos(tt);
            nodeY(numNodes,1) = rr * sin(tt);
            
            if irad == 1
                ss = pi/2 * Ra / (nn-1);
                if ii==1 || ii==nn
                    ss = ss * 0.5;
                end
                pp = P0 * ss * [cos(tt),sin(tt)];
                
                loaded_nodes(end+1,1) = numNodes;
                loaded_force(1:2,end+1) = pp;
            end
            
            %
            if irad == 1
                if ii < nn
                    numEdges = numEdges + 1;
                    edgePair(numEdges,:) = [numNodes+1,numNodes];
                end
            end
            if irad == nrad1
                if ii < nn
                    numEdges = numEdges + 1;
                    edgePair(numEdges,:) = [numNodes,numNodes+1];
                end
            end
            if ii == 1
                if irad < nrad1
                    numEdges = numEdges + 1;
                    edgePair(numEdges,:) = [numNodes,numNodes+nn];
                end
            end
            if ii == nn
                if irad < nrad1
                    numEdges = numEdges + 1;
                    edgePair(numEdges,:) = [numNodes+nn,numNodes];
                end
            end
        end
    end
    
    % setup Dirichlet BC
    % constrained on x-axis
    fixed_nodes_y = find(abs(nodeY)<1.0e-6);
    tmp = fdpmNodeDof(fixed_nodes_y,'y');
    dispBCDofs = [dispBCDofs; tmp];
    dispBCVals = [dispBCVals; zeros(size(tmp))];
    % constrained on y-axis
    fixed_nodes_x = find(abs(nodeX)<1.0e-6);
    tmp = fdpmNodeDof(fixed_nodes_x,'x');
    dispBCDofs = [dispBCDofs; tmp];
    dispBCVals = [dispBCVals; zeros(size(tmp))];
    
    % setup Loading
    tracBCDofs = fdpmNodeDof(loaded_nodes,'xy');
    tracBCVals = reshape(loaded_force, [],1);
    
elseif prob_type == 2
    % axisymmetric
    % simply generate a rectangular block
    
    % particle number in height direction
    nlen = round(H/h0);
    
    nrad1 = nrad + 1;
    nlen1 = nlen + 1;
    
    % used to find most inner particles for apply pressure
    loaded_nodes = [];
    loaded_force = [];
    % used to find lower/upper particles to constraint on x-axis
    fixed_nodes_y = [];
    fixed_nodes_x = [];
    
    for irad = 1:nrad1
    for jlen = 1:nlen1
        if 1
            rr = Ra + (irad-1)*h0;
            zz = 0 + (jlen-1)*h0;
            
            numNodes = numNodes + 1;
            
            nodeX(numNodes,1) = rr;
            nodeY(numNodes,1) = zz;
            
            nodeVol(numNodes,1) = 2*pi*rr * h0^2;
            
            if irad == 1
                pp = P0*(2*pi*Ra)*h0;
                if jlen==1 || jlen==nlen1
                    pp = pp * 0.5;
                end
                
                loaded_nodes(end+1,1) = numNodes;
                loaded_force(1:2,end+1) = [pp,0];
            end
            
            if jlen==1 || jlen==nlen1
                fixed_nodes_y(end+1,1) = numNodes;
            end
        end
    end
    end
    
    % setup Dirichlet BC
    % constrained on x-axis
    tmp = fdpmNodeDof(fixed_nodes_y,'y');
    dispBCDofs = [dispBCDofs; tmp];
    dispBCVals = [dispBCVals; zeros(size(tmp))];
    
    % setup Loading
    tracBCDofs = fdpmNodeDof(loaded_nodes,'xy');
    tracBCVals = reshape(loaded_force, [],1);
    
else
    error('Unknown prob_type');
end



