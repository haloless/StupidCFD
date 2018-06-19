% testVonMisesDriverUniaxial
% uniaxial test

maxiter = 25;
% tolrel = 1.0e-8;
% tolabs = 1.0e-8;
tol = 1.0e-8;

straindof = [1];
stressdof = 1:6; stressdof(straindof) = [];


maxstrain = 0.01;
maxstep = 100;

de11step = [];
de11step(end+1:end+100) = maxstrain / 100;
de11step(end+1:end+200) = -maxstrain / 100;
de11step(end+1:end+100) = maxstrain / 100;
% de11step(end+1:end+100) = -maxstrain / 100;

nstep = numel(de11step);

e11 = 0;

data = [];

for step = 1:nstep
    
    epsEold = epsE;
    epbarold = epbar;
    sigmayold = sigmay;
    
    e11old = e11;
    de11 = de11step(step);
    e11 = e11old + de11;
    
    % incremental strain
    % only the strain-controled directions are known
    % in this test, e11 is applied
    % while [e22,e33,e12,e23,e31] must be solved to achieve uniaxial stress state
    deps = zeros(6,1);
    deps(straindof) = de11;
    
    % newton iteration to achieve stress state
    conv = 0;
    for iter = 1:maxiter
        
        epsEtr = epsEold + deps;
        [dgam, sigma, epsE, epbar, sigmay, Dalg] = material3dVonMisesLinHard(epsEtr, epbarold, E,nu,sigmayold,H);
        
        rhs = zeros(6,1) - sigma;
        rhs(straindof) = 0;
        
        rnorm = norm(rhs);
        % disp(['iter=',int2str(iter),',rnorm=',num2str(rnorm)]);
        % check convergence
        if rnorm <= tol
            conv = 1;
            break;
        end
        
        ddeps = Dalg(stressdof,stressdof) \ rhs(stressdof);
        deps(stressdof) = deps(stressdof) + ddeps;
    end
    
    if ~conv
        error('Newton failed');
    end
    
    data(:,step+1) = [ e11, sigma(1) ];
    
end

figure;
plot(data(1,:),data(2,:),'x-');
xlabel('e11');
ylabel('sigma11');



