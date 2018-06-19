
//
R0 = 6.413e-3;
Rmid = R0 * 0.982;
L0 = 53.334e-3;
Lmid = L0 * 0.5;

//
dh1 = 0.4e-3;
dh2 = 2e-3;


Point(1) = {0,0,0, dh1};
Point(2) = {Rmid,0,0, dh1};
Point(3) = {R0,Lmid,0, dh2};
Point(4) = {0,Lmid,0, dh2};




//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line Loop(5) = {4, 1, 2, 3};
//+
Plane Surface(6) = {5};
//+
Physical Line("left") = {4};
//+
Physical Line("right") = {2};
//+
Physical Line("bottom") = {1};
//+
Physical Line("top") = {3};
//+
Physical Surface("bar") = {6};
