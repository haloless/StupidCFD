function [d,react] = fdpmSolve(K,f, bcdof,bcval)
%fdpmSolve
% stiffness matrix K
% external - internal force f
% 
% We remove Dirichlet BC
%

% free dof
nall = size(f,1);
freedof = (1:nall)';
% remove Dirichlet BC dof, this reduces the vector size
freedof(bcdof) = [];

% rhs
if isempty(bcval)
    % homogeneous BC
    rhs = f(freedof);
else
    rhs = f(freedof) - K(freedof,bcdof)*bcval;
end

%
sol = K(freedof,freedof) \ rhs;

% displacement as solution
d = zeros(size(f));
% fill solved part
d(freedof) = sol;
if ~isempty(bcval)
    % fill Dirichlet BC value
    d(bcdof) = bcval;
end

if nargout > 1
    % react force
    % this should only be nonzero for Dirichlet points
    react = K * d - f;
end

return
end

