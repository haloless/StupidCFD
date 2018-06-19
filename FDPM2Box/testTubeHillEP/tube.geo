
R1 = 1.0;
R2 = 2.0;

dh1 = 0.02;
dh2 = 0.2;


Point(1) = {0, 0, 0, dh1};
Point(2) = {R1, 0, 0, dh1};
Point(3) = {R2, 0, 0, dh2};
Point(4) = {0, R2, 0, dh2};
Point(5) = {0, R1, 0, dh1};

Line(1) = {2, 3};
Circle(2) = {3, 1, 4};
Line(3) = {4, 5};
Circle(4) = {5, 1, 2};

Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};

Physical Line("bottom") = {1};
Physical Line("left") = {3};
Physical Line("inner") = {4};
Physical Line("outer") = {2};
Physical Surface("tube") = {6};



