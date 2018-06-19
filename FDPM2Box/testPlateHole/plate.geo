
R = 1;
L = 5;


size1 = 0.05;
size2 = 0.5;


Point(1) = {0, 0, 0, size1};
Point(2) = {R, 0, 0, size1};
Point(3) = {0, R, 0, size1};
Point(4) = {L, 0, 0, size2};
Point(5) = {L, L, 0, size2};
Point(6) = {0, L, 0, size2};

Line(1) = {2, 4};
Line(2) = {4, 5};
Line(3) = {5, 6};
Line(4) = {6, 3};
Circle(5) = {3, 1, 2};
Line Loop(7) = {4, 5, 1, 2, 3};

Plane Surface(7) = {7};

//+
Physical Line("bottom") = {1};
//+
Physical Line("left") = {4};
//+
Physical Line("right") = {2};
//+
Physical Line("top") = {3};
//+
Physical Line("hole") = {5};
//+
Physical Surface("plate") = {7};
