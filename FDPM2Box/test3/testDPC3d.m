
clear;

par = struct();

% elastic
par.E = 100;
par.nu = 0.2;

E = par.E;
nu = par.nu;
De = E/(1+nu)/(1-2*nu) * [...
1-nu, nu, nu, 0, 0, 0;
nu, 1-nu, nu, 0, 0, 0;
nu, nu, 1-nu, 0, 0, 0;
0, 0, 0, 0.5-nu, 0, 0;
0, 0, 0, 0, 0.5-nu, 0;
0, 0, 0, 0, 0, 0.5-nu;
];

% DP
par.c0 = 0.1;

friction_angle = 20 / 180 * pi;
dilation_angle = 10 / 180 * pi;
% outer approx.
% par.eta    = 6/sqrt(3) * sin(friction_angle) / (3-sin(friction_angle));
% par.etabar = 6/sqrt(3) * sin(dilation_angle) / (3-sin(dilation_angle));
% par.xi     = 6/sqrt(3) * cos(friction_angle) / (3-sin(friction_angle));
% plain-strain approx.
par.eta    = 3 * tan(friction_angle) / sqrt(9 + 12*tan(friction_angle)^2);
par.etabar = 3 * tan(dilation_angle) / sqrt(9 + 12*tan(dilation_angle)^2);
par.xi     = 3 / sqrt(9 + 12*tan(friction_angle)^2);

% Cap
par.ptens = par.xi / par.eta * par.c0;
par.M = par.eta * sqrt(3);
par.beta = 0.25;
par.Hcap = 10;

if 0
    % set initial cap position at zero pressure
    % i.e. not compacted
    par.a0 = par.ptens / (1+par.beta);
else
    % assign inital cap position
    p0 = -0.1;
    par.a0 = (par.ptens-p0) / (1+par.beta);
end

% 
epsEP = zeros(6,1);
epsE = zeros(6,1);
epbar = 0;
alpha = 0;
sigma = zeros(6,1);
a = par.a0;

if 1
    hfig = figure;
    
    % plot DP envelop
    pp = -0.5:0.01:0.3;
    qq = -par.eta*pp + par.xi*par.c0;
    plot(pp,qq,'g-', pp,zeros(size(pp)),'k-', zeros(size(qq)),qq,'k-');
    xlabel('p'); ylabel('sqrt(J2)'); % ylabel('$$\sqrt{J_2}$$', 'Interpreter','latex');
    
    hold on;
    
    % plot cap
    plotCap(par, a, 'r-');
    
    % plot initial stress
    [p,s] = voigt3dPressShear(sigma);
    sj2 = sqrt(voigt3dJ2(s));
    plot(p,sj2,'x');
    
    hold off;
    % return;
end


% isotropic consolidation
if 1
    % testDPCIsoConsolidate;
end

% uniaxial strain
if 1
    testDPC3dUniaxialStrain;
end

% triaxial test
if 1
    % testDPCTriaxialShear;
end





