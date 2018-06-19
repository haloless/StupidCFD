% Build particle connectivity
% 
%

clear conn;

for i = 1:numNodes
    
    [neigh,rs,rx,ry,cutoff] = fdpmNeighborhood(i,nodePos,re);
    neigh2 = [neigh; i]; % particle included in the neighbourhood
    
    [N,Nx,Ny] = fdpmShapeMLS2(rx(neigh2),ry(neigh2),cutoff);
    
	% save connection
    conn(i).numNeigh = numel(neigh2);
	conn(i).neigh2 = neigh2'; % neighbor list is a row vector
    % conn(i).W = [w(:)', 0]; % weights, the last self weight is zero
    conn(i).N = N;
	conn(i).NX = Nx;
	conn(i).NY = Ny;
end







