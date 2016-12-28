
TwoPi = pi * 2;
ThreePi = pi * 3;
FourPi = pi * 4;
SixPi = pi * 6;
EightPi = pi * 8;


Delta = eye(3);

Epsil = zeros(3,3,3);
Epsil(1,2,3) = 1;
Epsil(2,3,1) = 1;
Epsil(3,1,2) = 1;
Epsil(3,2,1) = -1;
Epsil(2,1,3) = -1;
Epsil(1,3,2) = -1;
