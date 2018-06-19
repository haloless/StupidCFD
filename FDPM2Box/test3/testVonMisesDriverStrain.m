% testVonMisesDriverStrain
% strain controled

% maxstrain = 0.01;
maxstrain = 0.001;
maxstep = 100;

de11step = [];
de11step(end+1:end+100) = maxstrain / 100;
% de11step(end+1:end+200) = -maxstrain / 100;
% de11step(end+1:end+100) = maxstrain / 100;
% de11step(end+1:end+100) = -maxstrain / 100;

nstep = numel(de11step);

e11 = 0;

data = [];
for step = 1:nstep
    
    e11old = e11;
    epsEold = epsE;
    epbarold = epbar;
    sigmayold = sigmay;
    
    
    % e11 = maxstrain / maxstep * step;
    % de11 = maxstrain / maxstep;
    de11 = de11step(step);
    e11 = e11 + de11;
    
    deps = [ de11, -nu*de11, -nu*de11, 0, 0, 0 ]';
    
    epsEtr = epsE + deps;
    [dgam, sigma, epsE, epbar, sigmay, Dalg] = material3dVonMisesLinHard(epsEtr, epbar, E,nu,sigmay,H);
    
    % [~,s] = voigt3dVolDev(sigma);
    [p,s] = voigt3dPressShear(sigma);
    q = sqrt(3/2) * voigt3dNorm(s);
    data(:,step+1) = [ e11, q*sign(s(1)) ];
    
    if step == nstep
        % use finite difference to check Consistent Tangent Operator
        % Dalg = d(sigma_{n+1}) / d(eps_{trial})
        Ddiff = zeros(6);
        dd = 1.0e-6
        
        for dir = 1:6
            epsEtrp = epsEtr; epsEtrp(dir) = epsEtrp(dir) + dd;
            [~,sigmap,~,~,~,~] = material3dVonMisesLinHard(epsEtrp, epbarold, E,nu,sigmayold,H);
            epsEtrm = epsEtr; epsEtrm(dir) = epsEtrm(dir) - dd;
            [~,sigmam,~,~,~,~] = material3dVonMisesLinHard(epsEtrm, epbarold, E,nu,sigmayold,H);
            Ddiff(:,dir) = (sigmap-sigmam)./(dd*2);
        end
        
        Ddiff
        Dalg
    end
end

figure;
plot(data(1,:),data(2,:),'x-');
xlabel('e11');
ylabel('q');



