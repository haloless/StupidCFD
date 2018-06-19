%


% epsxmax = -0.01;
% epsymax = 0.0;
% epszmax = 0.0;

% epsxmax = -0.01;
% epsymax = 0.005;
% epszmax = 0.0;

% epsxmax = 0.01;
% epsymax = 0.01;
% epszmax = 0.01;

% epsxmax = 0.01;
% epsymax = -0.006;
% epszmax = -0.006;

% epsxmax = 0.0;
% epsymax = -0.01;
% epszmax = 0.0;

epsxmax = 0.006;
epsymax = -0.01;
epszmax = 0.006;



data = [];

nstep = 40;
for istep = 1:nstep
    
    epsEold = epsE;
    epbarold = epbar;
    alphaold = alpha;
    
    depsx = epsxmax / nstep;
    depsy = epsymax / nstep;
    depsz = epszmax / nstep;
    if prob_type == 1
        deps = [ depsx; depsy; 0; 0 ];
    elseif prob_type == 2
        deps = [ depsx; depsy; 0; depsz ];
    end
    
    epsEtr = epsE + deps;
    
    if 1
        figure(hfig);
        hold on;
        sigmatr = De * epsEtr;
        [ptr,str] = voigtPressShear(sigmatr);
        plot(ptr,voigtSqrtJ2(str),'.');
        hold off;
        pause;
    end
    
    [dgama,dgamb,sigma,epsE,epbar,alpha, Dalg,mtype] = materialDPC(par,prob_type, epsEtr,epbarold,alphaold);
    
    epsEP = epsEP + deps;
    % sigma
    
    if 0
        figure(hfig);
        hold on;
        testDPC3dReplotCurr;
        hold off;
        pause;
    end
    
    if 0
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




