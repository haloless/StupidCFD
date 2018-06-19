% Generate particles and determine loadings
%
%

if prob_type == 1
    % plane strain
    % generate a circular slice of the tube
    
    % used to find most inner particles for apply pressure
    loaded_nodes = [];
    loaded_force = [];
    
    for irad = 1:nrad
        if 1
            % 1/4 cylinder, constraint on x- and y-axis
            rr = Ra + (irad-0.5)*h0;
            ll = pi/2 * rr;
            nn = round(ll/h0)+1;
            vv = pi/4 * ((rr+0.5*h0)^2-(rr-0.5*h0)^2) / (nn-1);
            
            for ii = 1:nn
                tt = pi/2 / (nn-1) * (ii-1);
                
                numNodes = numNodes + 1;
                nodeX(numNodes,1) = rr * cos(tt);
                nodeY(numNodes,1) = rr * sin(tt);
                if ii==1 || ii==nn
                    nodeVol(numNodes,1) = vv;
                else
                    nodeVol(numNodes,1) = vv;
                end
                
                
                if irad == 1
                    pp = P0*(pi/2*Ra)/(nn-1) * [cos(tt),sin(tt)];
                    if ii==1 || ii==nn
                        pp = pp * 0.5;
                    end
                    
                    loaded_nodes(end+1,1) = numNodes;
                    loaded_force(1:2,end+1) = pp;
                end
            end
        else
            % full cylinder, choose fixed points lie on x- and y-axis 
            rr = Ra + (irad-0.5)*h0;
            ll = pi*2 * rr;
            nn = round(ll/h0);
            vv = pi * ((rr+0.5*h0)^2-(rr-0.5*h0)^2) / (nn);
            
            for ii = 1:nn
                tt = pi*2 / (nn) * (ii-1);
                
                numNodes = numNodes + 1;
                nodeX(numNodes,1) = rr * cos(tt);
                nodeY(numNodes,1) = rr * sin(tt);
                nodeVol(numNodes,1) = vv;
                
                if irad == 1
                    loaded_nodes(end+1,1) = numNodes;
                    loaded_force(1:2,end+1) = P0*(pi*2*Ra)/(nn) * [cos(tt),sin(tt)];
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
    
    
    % used to find most inner particles for apply pressure
    loaded_nodes = [];
    loaded_force = [];
    % used to find lower/upper particles to constraint on x-axis
    fixed_nodes_y = [];
    fixed_nodes_x = [];
    
    for irad = 1:nrad
    for jlen = 1:nlen
        if 1
            rr = Ra + (irad-0.5)*h0;
            zz = 0 + (jlen-0.5)*h0;
            
            numNodes = numNodes + 1;
            
            nodeX(numNodes,1) = rr;
            nodeY(numNodes,1) = zz;
            
            nodeVol(numNodes,1) = 2*pi*rr * h0^2;
            
            if irad == 1
                pp = P0*(2*pi*Ra)*h0;
                if jlen==1 || jlen==nlen
                    % pp = pp * 0.5;
                end
                
                loaded_nodes(end+1,1) = numNodes;
                loaded_force(1:2,end+1) = [pp,0];
            end
            
            if jlen==1 || jlen==nlen
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



