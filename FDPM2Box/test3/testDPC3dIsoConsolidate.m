%testDPCIsoConsolidate


presmax = -1;
% presmax = par.ptens;

tol = abs(presmax) * 1.0e-9;


data = [];

nstep = 20;
for istep = 1:nstep
    
    epsEold = epsE;
    epbarold = epbar;
    alphaold = alpha;
    
    % external load stress
    load = zeros(6,1);
    load(1:3) = presmax / nstep * istep;
    
    % strain increment
    deps = zeros(6,1);
    % deps(1:3) = 1.0e-6 * sign(presmax);
    
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
    
    if 1
        figure(hfig);
        hold on;
        testDPCReplotCurr;
        hold off;
        pause;
    end
end



