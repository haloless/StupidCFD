function [Dalgo,sigma,epsE,epflag] = materialVonMises(epsEtr, E,nu,fc, prob_type)
%materialVonMises: von Mises plastic, perfect (no hardening)
% INPUT
% * E, nu: Young and Poisson
% * epsEtr [11,22,12,33]
% * fc: if FC<0 return elastic state directly
%

assert(numel(epsEtr) == 4);

% set default values
if ~exist('prob_type','var')
    prob_type = 1;
end


tol = 1.0e-8;
maxiter = 25;

sqrt3 = sqrt(3);
bm1 = [ 1; 1; 0; 1 ];

if prob_type == 1
    ncomp = 3;
elseif prob_type == 2
    ncomp = 4;
end

% 4x4 elastic
De = E/(1+nu)/(1-2*nu) * [...
1-nu, nu,   0,      nu; 
nu,   1-nu, 0,      nu;
0,    0,    0.5-nu, 0; 
nu,   nu,   0,      1-nu];

% trial stress
sigma = De * epsEtr;
% set trial state
Dalgo = De(1:ncomp,1:ncomp);
epsE = epsEtr;
% flag
epflag = 0;

if fc < 0
    return % elastic, do not bother check plastic
end

% pressure
p = sum(sigma([1,2,4])) / 3;
% deviatoric
s = sigma - p.*bm1;
% J2 invariant
j2 = 0.5 * (s(1)^2 + s(2)^2 + s(4)^2 + 2*s(3)^2);

% yield
f = sqrt(3*j2) / fc - 1;

if (f > tol)
    % yield, elastoplastic
    epflag = 1;
    
    % residual vector
    b = zeros(5,1);
    b(end) = f;
    
    % plastic multiplier
    dgam = 0;
    
    ok = 0;
    for iter = 1:maxiter
        
        % d(J2)/d(sigma) = [s11,s22,2*s12,s33]
        dj2 = s;
        dj2(3) = 2*s(3);
        %
        ddj2 = zeros(4);
        ddj2([1,2,4],[1,2,4]) = eye(3) - ones(3)/3;
        ddj2(3,3) = 2;
        
        %
        df = sqrt3/2/fc / sqrt(j2) * dj2;
        ddf = sqrt3/2/fc* (-(dj2*dj2')/(2*j2^(3/2)) + ddj2/sqrt(j2));
        
        b = [ epsE-epsEtr+dgam*df; f ];
        
        if norm(b(1:4))<=tol && abs(b(5))<=tol
            ok = 1; break;
        end
        
        A = [ eye(4)+dgam*ddf*De, df; df'*De, 0 ];
        dx = -A \ b;
        
        % update strain and multiplier
        epsE = epsE + dx(1:4);
        dgam = dgam + dx(5);
        
        % update stress
        sigma = De * epsE;
        p = sum(sigma([1,2,4])) / 3;
        s = sigma - p.*bm1;
        j2 = 0.5 * (s(1)^2 + s(2)^2 + s(4)^2 + 2*s(3)^2);
        
        % check
        f = sqrt(3*j2) / fc - 1;
        
    end
    
    if ~ok
        error('material failed to converge');
    end
    
    % consistent tangent modulus
    Is = eye(4); 
    % Is(3,3) = 0.5; % do not set this
    P = (Is + dgam*De*ddf) \ De;
    Pdf = P * df;
    Dalgo = P - Pdf*Pdf'/(df'*Pdf);
    Dalgo = Dalgo(1:ncomp,1:ncomp);
end



return
end


