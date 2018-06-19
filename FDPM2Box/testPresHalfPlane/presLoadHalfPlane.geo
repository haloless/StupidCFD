
a = 1.0;
L = 3.0;

yprobe = 2.0;

dh1 = 0.2;
dh2 = 0.2;

Point(1) = {0, 0, 0, dh1};
Point(2) = {L, 0, 0, dh1};
Point(3) = {L, L, 0, dh2};
Point(4) = {0, L, 0, dh1};
Point(5) = {L-a, L, 0, dh2}; // 
Point(6) = {0, yprobe, 0, dh1}; // 
Point(7) = {L, yprobe, 0, dh1}; // 




Line(1) = {1, 2};
Line(2) = {2, 7};
Line(3) = {7, 3};
Line(4) = {3, 5};
Line(5) = {5, 4};
Line(6) = {4, 6};
Line(7) = {6, 1};
Line(8) = {6, 7};

Line Loop(9) = {7, 1, 2, -8};
Plane Surface(10) = {9};
Line Loop(11) = {6, 8, 3, 4, 5};
Plane Surface(12) = {11};


Physical Line("left") = {6, 7};
Physical Line("right") = {2, 3};
Physical Line("bottom") = {1};
Physical Line("topfree") = {5};
Physical Line("topload") = {4};
Physical Line("yline") = {8};
Physical Surface("plane") = {10, 12};





