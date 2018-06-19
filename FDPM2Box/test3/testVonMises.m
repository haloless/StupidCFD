
clear;


E = 210;
nu = 0.3;
sigmay0 = 0.24;
% sigmay0 = 10000;
% H = 0;
H = 10;

% elastic strain
epsE = zeros(6,1);
% yield stress
sigmay = sigmay0;
% accumulative plastic
epbar = 0;


% strain controled
testVonMisesDriverStrain;

% uniaxial
% testVonMisesDriverUniaxial;




