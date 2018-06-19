function [j2,dj2,ddj2] = voigt3dJ2(s)
%voigt3dJ2
% Input: s is the deviatoric part of stress sigma
% Output: J2, dJ2/dsigma, dd(J2)/dd(sigma)

% J2
% notice the factor 2
j2 = 0.5 * (sum(s(1:3).^2) + 2*sum(s(4:6).^2));

if nargout > 1 
    % dJ2/dsig
    dj2 = zeros(6,1);
    dj2(1:3) = s(1:3);
    dj2(4:6) = 2*s(4:6); % notice the factor 2
end

if nargout > 2
    % ddJ2/ddsig
    ddj2 = zeros(6);
    ddj2(1:3,1:3) = eye(3) - 1/3*ones(3);
    ddj2(4:6,4:6) = 2*eye(3); % notice the factor 2
end


return
end


