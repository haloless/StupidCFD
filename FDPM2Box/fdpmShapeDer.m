function [dN,dNX,dNY,dNXX,dNXY,dNYY] = fdpmShapeDer(rx,ry,wgt)
%fdpmShapeDer
% Calculate derivative shape function.
%
% Input:
% rx,ry: contains relative positions
% wgt: weights
%
% Output:
%


% support polynomial
P = fdpmShapeBasis(rx,ry);

% moment matrix
Amat = P' * diag(wgt) * P;
dN = (Amat' \ P') * diag(wgt);

% append the self part
dN = [ dN, -sum(dN,2) ];

%
if nargout > 1
	dNX  = dN(1,:);
	dNY  = dN(2,:);
	dNXX = dN(3,:);
	dNXY = dN(4,:);
	dNYY = dN(5,:);
end

% % save shape function
% dNX(i,neigh) = dN(1,:);
% dNY(i,neigh) = dN(2,:);
% dNX(i,i) = -sum(dNX(i,neigh));
% dNY(i,i) = -sum(dNY(i,neigh));

% dNXX(i,neigh) = dN(3,:);
% dNXX(i,i) = -sum(dNXX(i,neigh));
% dNXY(i,neigh) = dN(4,:);
% dNXY(i,i) = -sum(dNXY(i,neigh));
% dNYY(i,neigh) = dN(5,:);
% dNYY(i,i) = -sum(dNYY(i,neigh));



return
end

