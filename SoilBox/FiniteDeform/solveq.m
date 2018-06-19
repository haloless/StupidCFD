function [d,R] = solveq(K,f,bc)

fdof = [1:size(K)]';

d = zeros(size(fdof));

% dof with known BC
fbc = bc(:,1);

% remove known BC, this reduces the vector size 
fdof(bc(:,1)) = [];

rhs = f(fdof) - K(fdof,fbc)*bc(:,2);

sol = K(fdof,fdof) \ rhs;

d(fbc) = bc(:,2);
d(fdof) = sol;

% residual
R = K * d - f;

return
end


