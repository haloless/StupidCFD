function [j2,dj2,ddj2] = voigt3dJ2(s)
%voigtJ2
% Input: s is the deviatoric part of stress sigma
% Output: J2, dJ2/dsigma, dd(J2)/dd(sigma)

assert(numel(s) == 4);

% J2
% notice the factor 2
j2 = 0.5 * (s(1)^2 + s(2)^2 + s(4)^2 + 2*s(3)^2);

if nargout > 1 
    % dJ2/dsig
    dj2 = s;
    dj2(3) = s(3)*2; % notice the factor 2
end

if nargout > 2
    % ddJ2/ddsig
    ddj2 = zeros(4);
    ddj2([1,2,4],[1,2,4]) = eye(3) - 1/3*ones(3);
    ddj2(3,3) = 2; % notice the factor 2
end


return
end


