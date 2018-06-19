%testDPCTriaxialShear
% stress state [sigxx,sigyy,sigzz,0,0,0]
% strain state [epsxx,epsyy,epszz,0,0,0]
% take sigyy=sigzz, epsyy=epszz
%
% x-axis is strain controlled
% y- and z-axis are stress controlled
%
% 1. isotropic consolidate to given pressure
% 2. compact in x, while keep stress in y&z


%
pcons = -0.1;

%
% exmax = -0.06;
exmax = -0.2;


tol = abs(pcons) * 1.0e-9;

%
% 1. isotropic consolidation
% 
disp('Begin isotropic consolidation');

nstep = 20;
for istep = 1:nstep
    
    epsEold = epsE;
    epbarold = epbar;
    alphaold = alpha;
    
    % external load stress
    load = zeros(6,1);
    load(1:3) = pcons / nstep * istep;
    
    % strain increment
    deps = zeros(6,1);
    
    % newtons iteration to solve strain increment to achieve applied stress
    conv = 0;
    for iter = 1:25
        
        epsEtr = epsEold + deps;
        [dgama,dgamb,sigma,epsE,epbar,alpha, Dalg,mtype] = material3dDPCap(par, epsEtr,epbarold,alphaold);
        
        rhs = load - sigma;
        
        rnorm = norm(rhs);
        if rnorm <= tol
            conv = 1; break;
        end
        
        ddeps = Dalg \ rhs;
        deps = deps + ddeps;
    end
    
    if ~conv
        error('Newton failed');
    end
    
    % accumulate full strain for use
    epsEP = epsEP + deps;
    
    if 1
        figure(hfig);
        hold on;
        testDPCReplotCurr;
        hold off;
        pause;
    end
end

disp('End isotropic consolidation');
% return

%
% 2. triaxial test
% 
disp('Begin triaxial test');

nstep = 100;

if abs(exmax) <= abs(epsEP(1))
    error('too much consolidation');
else
    depsxx = (exmax - epsEP(1)) / nstep;
    disp(['epsxx=',num2str(epsEP(1)),'; epsxxmax=',num2str(exmax)]);
    disp(['depsxx=',num2str(depsxx)]);
end

for istep = 1:nstep
    
    epsEold = epsE;
    epbarold = epbar;
    alphaold = alpha;
    
    % external load stress
    load = zeros(6,1);
    load(2:3) = pcons;
    
    % strain increment
    deps = zeros(6,1);
    deps(1) = depsxx;
    
    % newtons iteration to solve strain increment to achieve applied stress
    conv = 0;
    for iter = 1:25
        
        epsEtr = epsEold + deps;
        [dgama,dgamb,sigma,epsE,epbar,alpha, Dalg,mtype] = material3dDPCap(par, epsEtr,epbarold,alphaold);
        
        rhs = load - sigma;
        rhs(1) = 0;
        
        rnorm = norm(rhs);
        if rnorm <= tol
            conv = 1; break;
        end
        
        ddeps = Dalg(2:6,2:6) \ rhs(2:6);
        deps(2:6) = deps(2:6) + ddeps;
    end
    
    if ~conv
        error('Newton failed');
    end
    
    % accumulate full strain for use
    epsEP = epsEP + deps;
    
    if 1
        figure(hfig);
        hold on;
        testDPCReplotCurr;
        hold off;
        pause;
    end
end

disp('End triaxial test');

