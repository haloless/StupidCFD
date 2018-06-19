function [dgama,dgamb,sigma,epsE,epbar,alpha,Dalgo,mtype] = materialDPC(par,probtype, epsEtr, epbarold,alphaold)
%materialDPC: 2D version (plain strain / axisymmetric)
% Standard Drucker-Prager Cap model
% DP envelop: perfect plastic (cohesion c = const)
% Cap: hardening
% 
% MTYPE=0, elastic
% MTYPE=1, cone
% MTYPE=2, apex
% MTYPE=3, cap
% MTYPE=4, cone-cap intersection
%

% check input trial elastic strain
assert(numel(epsEtr) == 4);

if probtype == 1 % plane stain
    ncomp = 3;
elseif probtype == 2 % axisymmetric
    ncomp = 4;
else
    error('Unsupported probtype = %d', probtype);
end


bm1 = [ 1; 1; 0; 1 ];

sqrt2 = sqrt(2.0);
sqrt3 = sqrt(3.0);

tol = 1.0e-9;

% parameters
% young and poisson
E = par.E;
nu = par.nu;

% DP parameters
eta = par.eta;
etabar = par.etabar;
xi = par.xi;
c0 = par.c0;

% Cap parameters
ptens = par.ptens;
beta = par.beta;
M = par.M;

% elastic modulus
K = E / 3 / (1-2*nu);
G = E / 2 / (1+nu);

De = E/(1+nu)/(1-2*nu) * [...
1-nu, nu,   0,      nu; 
nu,   1-nu, 0,      nu;
0,    0,    0.5-nu, 0; 
nu,   nu,   0,      1-nu];

% elastic trial state
epsE = epsEtr;
sigma = De * epsEtr;
dgama = 0;
dgamb = 0;
epbar = epbarold;
alpha = alphaold;

% return-mapping type
mtype = -1;


%
[ptr,str] = voigtPressShear(sigma);
% J2 invariant
j2tr = voigtJ2(str);
% sqrt(J2)
sj2tr = sqrt(j2tr);
% q(s)
qtr = sqrt3 * sj2tr;


% check yield
% DP
fa = (sj2tr + eta*ptr) / (xi*c0) - 1;

% cap
aold = apresCap(par, alphaold);
fb = sqrt((ptr-ptens+aold)^2/beta^2 + qtr^2/M^2) / aold - 1;


if (ptr>=ptens-aold && fa<=tol) || (ptr<ptens-aold && fb<=tol) || true
    % elastic, do nothing
    mtype = 0;
    
else % if fa>tol || (fb>tol && ptr<ptens-aold) % plastic
    
    if fa > tol
        % 1. try return to Cone
        mtype = 1;
        [dgama,p,s,epbar] = retmapDPCone(par, ptr,str, epbarold, K,G,eta,etabar,xi,c0);
        
        % check apex 
        if sj2tr - dgama*G < 0
            % 2. return to Apex
            mtype = 2;
            [dgama ,p,s,epbar] = retmapDPApex(par, ptr,str, epbarold, K,G,eta,etabar,xi,c0);
            
        else
            % check cap
            q = sqrt(3 * voigtJ2(s));
            fbnew = sqrt((p-ptens+aold)^2/beta^2 + q^2/M^2) / aold - 1;
            if fbnew>tol && p<ptens-aold
                % invalid map, will try cap
                % reset plastic variables
                mtype = -1;
                
                dgama = 0;
                epbar = epbarold;
            end
        end
    end
    
    if mtype == -1
        assert(fb>=0);
        
        % 3. try return to Cap
        mtype = 3;
        [dgamb,p,s,alpha] = retmapCap(par, ptr,str, alphaold, K,G,ptens,beta,M);
        
        % check cap validity
        a = apresCap(par, alpha);
        if p > ptens - a
            % 4. return to cone-cap intersection
            mtype = 4;
            % error('Cone/Cap not implemented');
            [dgama,dgamb,p,s,epbar,alpha] = retmapCorner(par, ptr,str,epbarold,alphaold, K,G, eta,etabar,xi,c0, ptens,beta,M);
        end
    end
    
    assert(mtype ~= -1);
    
    % update stress
    sigma = s + p.*bm1;
    % elastic strain
    epsE = De \ sigma;
    
end

%
% compute Consistent Tangent
%

% hardening rate
% NOTE we assume perfect plasticity here
Hdp = 0;

if mtype == 0 % elastic
    
    Dalgo = De;
    
elseif mtype == 1 % cone
    dgam = dgama;
    
    % volumetric
    epsEtrv = sum(epsEtr([1,2,4]));
    
    % deviatoric
    epsEtrd = epsEtr - (epsEtrv/3).*bm1;
    % recover the [2*e12]
    epsEtrd(3) = 0.5 * epsEtrd(3);
    
    % norm
    epsEtrdnorm = sqrt(epsEtrd(1)^2 + epsEtrd(2)^2 + epsEtrd(4)^2 + 2*epsEtrd(3)^2);
    
    if epsEtrdnorm > 0
        ninv = 1 / epsEtrdnorm;
    else
        ninv = 0;
    end
    
    % edev / |edev|
    unidev = epsEtrd .* ninv;
    
    aux = 1.0 / (G + K*eta*etabar + xi^2*Hdp);
    
    afact = 2*G * (1 - dgam/sqrt(2)*ninv);
    bfact = 2*G * (dgam/sqrt(2)*ninv - G*aux);
    cfact = -sqrt(2) * G * aux * K;
    dfact = K * (1 - K*eta*etabar*aux);
    
    aaa = afact .* diag([1, 1, 1/2, 1]);
    bbb = bfact .* unidev * unidev';
    ccc = cfact .* (eta.*unidev*bm1' + etabar.*bm1*unidev');
    ddd = (dfact - afact/3) .* bm1*bm1';
    
    Dalgo = aaa + bbb + ccc + ddd;
    
elseif mtype == 2 % apex
    dgam = dgama;
    acoef = xi / etabar;
    bcoef = xi / eta;
    
    % perfect plasticity, this is essentially zero
    factor = K * (1 - K / (K+acoef*bcoef*Hdp));
    
    Dalgo = factor .* diag(bm1);
elseif mtype == 3 % cap
    dgam = dgamb;
    
    [a,aH] = apresCap(par, alpha);
    
    [~,dj2,ddj2] = voigtJ2(s);
    Ndev = 3/M^2 * dj2;
    Nvol = 2/beta^2 * (p-ptens+a);
    Nvec = Ndev + (Nvol/3).*bm1;
    dNddsigma = 3/M^2 * ddj2;
    dNvdsigma = 2/beta^2 / 3 * bm1;
    dNdsigma = dNddsigma + 1/3*(bm1*dNvdsigma');
    dNdalpha = 2/beta^2 / 3 * aH .* bm1;
    dfdalpha = (Nvol - 2*a) * aH;
    
    Cmat = [ inv(De)+dgam*dNdsigma, dgam*dNdalpha, Nvec;
    dgam*dNvdsigma', 1 + 2/beta^2*aH*dgam, Nvol;
    (Nvec)', dfdalpha, 0];
    
    Dmat = inv(Cmat);
    
    Dalgo = Dmat(1:4,1:4);
elseif mtype == 4
    
    b = beta;
    
    [a,aH] = apresCap(par, alpha);
    
    [~,dj2,ddj2] = voigtJ2(s);
    [~,dsj,ddsj] = voigtSqrtJ2(s);
    
    Na = dsj + (1/3*etabar).*bm1;
    dNadsigma = ddsj;
    
    Nb = 3/M^2*dj2;
    % dNbdsigma = 3/M^2*ddj2 + 2/b^2/9*(bm1*bm1');
    % dNbdalpha = 2/b^2 * aH /3 .* bm1;
    dNbdsigma = 3/M^2*ddj2;
    dNbdalpha = zeros(4,1);
    
    dalpha = [ zeros(1,4), 1, etabar, 0 ];
    % dalpha = [ dgamb*2/b^2/3*bm1', (1+dgamb*2/b^2*aH), etabar, 0 ];
    
    dfadsigma = dsj + (1/3*eta).*bm1;
    
    dfbdalpha = -2*a*aH;
    
    Cmat = [ inv(De)+dgama*dNadsigma+dgamb*dNbdsigma, dgamb*dNbdalpha, Na, Nb;
    dalpha;
    dfadsigma', 0, 0, 0;
    Nb', dfbdalpha, 0, 0];
    
    Dmat = inv(Cmat);
    
    Dalgo = Dmat(1:4,1:4);
else
    error('Invalid MTYPE');
end

% crop Dalgo to correct size according to problem type
Dalgo = Dalgo(1:ncomp,1:ncomp);


return
end

function [a,da] = apresCap(par, alpha)
%apresCap 
% hardening force a = a(alpha)
% derivative da/dalpha

azero = par.a0;
Hcap = par.Hcap;

% linear equation
a = azero + Hcap * alpha;
da = Hcap;

return
end



function [dgam,p,s,epbar] = retmapDPCone(par, ptr,str, epbarold, K,G,eta,etabar,xi,c0)
    
    % J2 invariant
    j2tr = voigtJ2(str);
    % sqrt(J2)
    sj2tr = sqrt(j2tr);
    
    % plastic multiplier
    dgam = (sj2tr + eta*ptr - xi*c0) / (G + eta*K*etabar);
    assert(dgam >= 0);
    
    % cumulative plastic strain
    depbar = xi * dgam;
    epbar = epbarold + depbar;
    
    % stress
    s = (1 - dgam*G/sj2tr) * str;
    p = ptr - dgam*K*etabar;
    
    return
end

function [dgam,p,s,epbar] = retmapDPApex(par, ptr,str, epbarold, K,G,eta,etabar,xi,c0)
    
    if eta==0 || etabar==0
        error('eta=0 or etabar=0 cannot return to apex');
    end
    
    acoef = xi / etabar;
    bcoef = xi / eta;
    
    % incremental volumetric plastic strain
    depv = (ptr - bcoef*c0) / K;
    
    % plastic multiplier
    dgam = depv / etabar;
    assert(dgam >= 0);
    
    % accumultive plastic strain
    depbar = acoef * depv;
    epbar = epbarold + depbar;
    
    % stress, only hydrostatic components
    p = ptr - K*depv;
    s = zeros(size(str));
    
    return
end

function [dgam,p,s,alpha] = retmapCap(par, ptr,str, alphaold, K,G,ptens,b,M)
    
    % aold = alphaold * Hslope + azero;
    % aold = apresCap(alphaold);
    qtr = sqrt(3 * voigtJ2(str));
    
    p = ptr;
    s = str;
    q = qtr;
    dgam = 0;
    alpha = alphaold;
    
    [a,aH] = apresCap(par,alpha);
    
    
    % solve (dgamma,alpha) by Newton
    tol = 1.0e-10;
    conv = 0;
    maxiter = 25;
    for iter = 1:maxiter
        pbar = p - ptens + a;
        
        rhs = zeros(2,1);
        rhs(1) = pbar^2/b^2 + q^2/M^2 - a^2;
        rhs(2) = alpha - alphaold + dgam*2/b^2*pbar;
        
        mat = zeros(2,2);
        mat(1,1) = -12*G*q^2/M^2/(M^2+6*G*dgam);
        mat(1,2) = 2/b^2*pbar*(K+aH) - 2*a*aH;
        mat(2,1) = 2/b^2*pbar;
        mat(2,2) = 1 + dgam*2/b^2*(K+aH);
        
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
        [a,aH] = apresCap(par,alpha);
        
    end
    
    if (~conv)
        error('Cam-Clay Newton iteration failed');
    end
    
    % check plastic multiplier
    assert(dgam >= 0);
    
    return
end

function [dgama,dgamb,p,s,epbar,alpha] = retmapCorner(par,ptr,str,epbarold,alphaold, K,G, eta,etabar,xi,c0, ptens,b,M)
    %
    
    sj2tr = sqrt(voigtJ2(str));
    qtr = sqrt(3) * sj2tr;
    
    dgama = 0;
    dgamb = 0;
    p = ptr;
    s = str;
    q = qtr;
    epbar = epbarold;
    alpha = alphaold;
    
    % solve (dot_gamma_a, dot_gamma_b) by Newton
    tol = 1.0e-10;
    conv = 0;
    maxiter = 25;
    for iter = 1:maxiter
        
        [a,aH] = apresCap(par, alpha);
        pbar = p - ptens;
        
        rhs = zeros(2,1);
        rhs(1) = q + sqrt(3)*eta*pbar;
        rhs(2) = (pbar+a)^2/b^2 + q^2/M^2 - a^2;
        
        sqrt3 = sqrt(3);
        M2 = M^2;
        b2 = b^2;
        gmb = G / (M2 + 6*G*dgamb);
        
        mat = zeros(2,2);
        mat(1,1) = -sqrt3*gmb*M2 - sqrt3*eta*K*etabar;
        mat(1,2) = -6*gmb*q;
        mat(2,1) = 2/b2*(pbar+a)*(-K*etabar-aH*etabar) - 2*sqrt3*gmb*q + 2*aH*etabar*a;
        mat(2,2) = -12 * gmb * q^2/M2;
        
        rnorm = norm(rhs);
        if rnorm <= tol
            conv = 1; break;
        end
        
        sol = mat \ (-rhs);
        dgama = dgama + sol(1);
        dgamb = dgamb + sol(2);
        
        p = ptr - K*etabar*dgama;
        s = M^2/(M^2+6*G*dgamb) * (1-G*dgama/sj2tr) * str;
        q = M^2/(M^2+6*G*dgamb) * (qtr - sqrt(3)*G*dgama);
        
        epbar = epbarold + xi*dgama;
        alpha = alphaold - etabar*dgama;
    end
    
    if (~conv)
        error('Corner newton failed');
    end
    
    assert(dgama >= 0);
    assert(dgamb >= 0);
    
    
    return
end





