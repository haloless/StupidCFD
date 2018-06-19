function [dgam,sigma,epsE,alpha,Dalg] = material3dCamClay(par, epsEtr,alphaold)
%material3dCamClay
% Assume linear hardening a(alpha) = a0 + H*alpha
% 


tol = 1.0e-12;
sqrt3 = sqrt(3.0);

bm1 = [ 1; 1; 1; 0; 0; 0 ];

bms = [ 1; 1; 1; 1/2; 1/2; 1/2 ];
Idev = diag([1, 1, 1, 2, 2, 2]) - (1/3).*(bm1*bm1');

%
E = par.E;
nu = par.nu;

%
ptens = par.ptens;
beta = par.beta;
M = par.M;
Hslope = par.Hslope;
azero = par.a0;

%
K = E / 3 / (1-2*nu);
G = E / 2 / (1+nu);

De = (K-2/3*G)/nu .* [...
1-nu, nu, nu, 0, 0, 0;
nu, 1-nu, nu, 0, 0, 0;
nu, nu, 1-nu, 0, 0, 0;
0, 0, 0, 1/2-nu, 0, 0;
0, 0, 0, 0, 1/2-nu, 0;
0, 0, 0, 0, 0, 1/2-nu;
];

%
% trial 
%

% epsEtrv = sum(epsEtr(1:3));
% epsEtrd = epsEtr - (epsEtrv/3).*bm1;
% ptr = K * epsEtrv;
% str = 2*G * bms .* epsEtrd;

sigmatr = De * epsEtr;
[ptr,str] = voigt3dPressShear(sigmatr);

% qtr = sqrt(3/2) * sqrt(sum(str(1:3).^2) + 2*sum(str(4:6).^2));
qtr = sqrt3 * voigt3dSqrtJ2(str);

% set elastic trial states
sigma = sigmatr;
epsE = epsEtr;
alpha = alphaold;
dgam = 0;
Dalg = De;

mtype = 0;

%
aold = alphaold * Hslope + azero;
if ptr >= ptens-aold
    b = 1;
    mtype = 1;
else
    b = beta;
    mtype = 2;
end

f = 1/b^2 * (ptr-ptens+aold)^2 + 1/M^2 * qtr^2 - aold^2;
if (f / aold^2 > tol)
    % plastic
    
    p = ptr;
    s = str;
    q = qtr;
    a = aold;
    
    conv = 0;
    maxiter = 25;
    for iter = 1:maxiter
        pbar = p - ptens + a;
        
        rhs = zeros(2,1);
        rhs(1) = pbar^2/b^2 + q^2/M^2 - a^2;
        rhs(2) = alpha - alphaold + dgam*2/b^2*pbar;
        
        mat = zeros(2);
        mat(1,1) = -12*G*q^2/M^2/(M^2+6*G*dgam);
        mat(1,2) = 2/b^2*pbar*(K+Hslope) - 2*a*Hslope;
        mat(2,1) = 2/b^2*pbar;
        mat(2,2) = 1 + dgam*2/b^2*(K+Hslope);
        
        rnorm = norm(rhs);
        if rnorm <= tol
            conv = 1; break;
        end
        
        sol = mat \ (-rhs);
        dgam = dgam + sol(1);
        alpha = alpha + sol(2);
        
        % update 
        s = M^2/(M^2+6*G*dgam) * str;
        q = M^2/(M^2+6*G*dgam) * qtr;
        p = ptr + K*(alpha-alphaold);
        a = alpha * Hslope + azero;
        
    end
    
    if (~conv)
        error('Cam-Clay Newton iteration failed');
    end
    
    % check plastic multiplier
    assert(dgam >= 0);
    
    % check yield surface
    if mtype == 1
        % assert(p >= ptens-a);
    elseif mtype == 2
        % assert(p < ptens-a, '(p=%f) < (ptens-a=%f)',p,ptens-a);
    end
    
    % stress
    sigma = s + p.*bm1;
    
    % strain
    % epsEd = s ./ (2*G*bms);
    % epsEv = p ./ K;
    % epsE = epsEd + (epsEv/3).*bm1;
    epsE = De \ sigma;
    
    %
    % compute cosistent tangent
    %
    
    Ce = inv(De);
    
    if 0
        % Ndev = 3/M^2 * s;
        Ndev = 3/M^2 * s ./ bms;
        Nvol = 2/b^2 * (p-ptens+a);
        Nvec = Ndev + (Nvol/3).*bm1;
        dNddsigma = 3/M^2 * Idev; 
        dNvdsigma = 2/(3*b^2) * bm1;
        dNdsigma = dNddsigma + 1/3*(bm1*dNvdsigma');
        dNdalpha = 2/(3*b^2) * Hslope .* bm1;
        dfdalpha = (Nvol - 2*a) * Hslope;
    else
        [~,dj2,ddj2] = voigt3dJ2(s);
        Ndev = 3/M^2 * dj2;
        Nvol = 2/b^2 * (p-ptens+a);
        Nvec = Ndev + (Nvol/3).*bm1;
        dNddsigma = 3/M^2 * ddj2;
        dNvdsigma = 2/b^2 / 3 * bm1;
        dNdsigma = dNddsigma + 1/3*(bm1*dNvdsigma');
        dNdalpha = 2/b^2 / 3 * Hslope .* bm1;
        dfdalpha = (Nvol - 2*a) * Hslope;
    end
    
    % Ndev
    % Nvol
    
    if 1
        % by full equations
        Cmat = [ Ce+dgam*dNdsigma, dgam*dNdalpha, Nvec;
        dgam*dNvdsigma', 1 + 2/b^2*Hslope*dgam, Nvol;
        (Nvec)', dfdalpha, 0];
        Dmat = inv(Cmat);
        Dalg = Dmat(1:6,1:6);
    else
        %
        Asub = Ce + dgam*dNdsigma;
        Bsub = [dgam*dNdalpha, Nvec];
        Csub = [dgam*dNvdsigma'; (Nvec)'];
        Dsub = [1 + 2/b^2*Hslope*dgam, Nvol; dfdalpha, 0];
        
        Ainv = inv(Asub);
        Dinv = inv(Dsub - Csub*Ainv*Bsub);
        Dalg = Ainv + Ainv*Bsub*Dinv*Csub*Ainv;
    end
end




return
end


