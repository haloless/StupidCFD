%

% epsxmax = -0.01;
% epsymax = -0.01;
% epszmax = -0.01;


epsxmax = -0.01;
epsymax = 0.0;
epszmax = 0.0;

% epsxmax = -0.01;
% epsymax = 0.005;
% epszmax = 0.0;

% epsxmax = 0.02;
% epsymax = 0.01;
% epszmax = 0.01;


% epsxmax = 0.01;
% epsymax = -0.006;
% epszmax = -0.006;

% epsxmax = 0.0;
% epsymax = -0.01;
% epszmax = 0.0;

% epsxmax = -0.01;
% epsymax = 0.006;
% epszmax = 0.006;


data = [];

nstep = 40;
for istep = 1:nstep
    
    epsEold = epsE;
    epbarold = epbar;
    alphaold = alpha;
    
    deps = [ epsxmax; epsymax; epszmax; 0; 0; 0 ] ./ nstep;
    
    epsEtr = epsE + deps;
    
    if 1
        figure(hfig);
        hold on;
        sigmatr = De * epsEtr;
        [ptr,str] = voigt3dPressShear(sigmatr);
        plot(ptr,voigt3dSqrtJ2(str),'.');
        hold off;
        pause;
    end
    
    [dgama,dgamb,sigma,epsE,epbar,alpha, Dalg,mtype] = material3dDPCap(par, epsEtr,epbarold,alphaold);
    
    epsEP = epsEP + deps;
    % sigma
    
    if 1
        figure(hfig);
        hold on;
        testDPC3dReplotCurr;
        hold off;
        pause;
    end
    
    if 1
        % use finite difference to check Consistent Tangent Operator
        % Dalg = d(sigma_{n+1}) / d(eps_{trial})
        Ddiff = zeros(6);
        dd = 1.0e-6;
        
        for dir = 1:6
            
            epsEtrp = epsEtr; epsEtrp(dir) = epsEtrp(dir) + dd;
            [~,~,sigmap] = material3dDPCap(par, epsEtrp,epbarold,alphaold);
            
            epsEtrm = epsEtr; epsEtrm(dir) = epsEtrm(dir) - dd;
            [~,~,sigmam] = material3dDPCap(par, epsEtrm,epbarold,alphaold);
            
            Ddiff(:,dir) = (sigmap-sigmam)./(dd*2);
        end
        
        Ddiff
        Dalg
        
        pause;
    end
end




