function [dgam,sigma,epsE,epbar,sigmay,Dalg] = material3dVonMisesLinHard(epsEtr, epbar, E,nu,sigmay,H)
%material3dVonMisesLinHard
% 

tol = 1.0e-9;

bm1 = [1;1;1;0;0;0];

% bulk
K = E / 3 / (1-2*nu);
G = E / 2 / (1+nu);

% elastic
De = E/(1+nu)/(1-2*nu) * [...
1-nu, nu, nu, 0, 0, 0;
nu, 1-nu, nu, 0, 0, 0;
nu, nu, 1-nu, 0, 0, 0;
0, 0, 0, 0.5-nu, 0, 0;
0, 0, 0, 0, 0.5-nu, 0;
0, 0, 0, 0, 0, 0.5-nu;
];

% trial state
epsE = epsEtr;
% trial stress
sigma = De * epsEtr;
% no plastic
dgam = 0;
% elastic tangent
Dalg = De;

% take deviatoric
% p = mean(sigma(1:3));
% s = sigma - p.*bm1;
[p,s] = voigt3dPressShear(sigma);
% J2 invariant
% j2 = 0.5 * voigt3dNorm(s)^2;
j2 = voigt3dJ2(s);
% q
q = sqrt(3*j2);

% check yield
f = q/sigmay - 1;

if (f > tol)
    % plastic
    
    %
    % return mapping 
    %
    dgam = (q - sigmay) / (3*G+H);
    assert(dgam >= 0);
    
    epbar = epbar + dgam;
    
    % update stress
    p = p;
    s = (1 - dgam*3*G/q) .* s;
    sigma = s + p.*bm1;
    
    epsE = De \ sigma;
    
    sigmay = sigmay + H*dgam;
    
    %
    % calculate consistent tangent
    %
    
    if 0
        % J2
        j2 = 0.5 * voigt3dNorm(s)^2;
        % dJ2/dsig
        dj2 = zeros(6,1);
        dj2(1:3) = s(1:3); dj2(4:6) = 2*s(4:6);
        % ddJ2/ddsig
        ddj2 = zeros(6);
        ddj2(1:3,1:3) = eye(3) - 1/3*ones(3);
        ddj2(4:6,4:6) = 2*eye(3);
        
        % dF/dsig
        df = sqrt(3)/2 / sqrt(j2) * dj2;
        % ddF/ddsig
        ddf = sqrt(3)/2 * (-(dj2*dj2')/(2*j2^(3/2)) + ddj2/sqrt(j2));
    else
        [sj2,dsj2,ddsj2] = voigt3dSqrtJ2(s);
        df = sqrt(3) * dsj2;
        ddf = sqrt(3) * ddsj2;
    end
    
    P = (eye(6) + dgam*De*ddf) \ De;
    Dalg = P - P*df*df'*P' / (df'*P*df+H);
end


return
end



