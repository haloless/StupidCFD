
L = 10;
C = 1;
D = C * 2;

dh1 = 0.2;
dh2 = 0.2;
dh3 = 0.2;


Point(1) = {0, -C, 0, dh1};
Point(2) = {L, -C, 0, dh1};
Point(3) = {L, C, 0, dh2};
Point(4) = {0, C, 0, dh2};
Point(5) = {L/2, -C, 0, dh3};
Point(6) = {L/2, C, 0, dh3};
Point(7) = {0, 0, 0, dh1};
Point(8) = {L, 0, 0, dh2};
Point(9) = {L/2, 0, 0, dh3};




//+
Line(1) = {1, 5};
//+
Line(2) = {5, 2};
//+
Line(3) = {2, 8};
//+
Line(4) = {8, 3};
//+
Line(5) = {3, 6};
//+
Line(6) = {6, 4};
//+
Line(7) = {4, 7};
//+
Line(8) = {7, 1};
//+
Line(9) = {7, 9};
//+
Line(10) = {9, 8};
//+
Line(11) = {5, 9};
//+
Line(12) = {9, 6};
//+
Line Loop(13) = {1, 11, -9, 8};
//+
Plane Surface(14) = {13};
//+
Line Loop(15) = {2, 3, -10, -11};
//+
Plane Surface(16) = {15};
//+
Line Loop(17) = {10, 4, 5, -12};
//+
Plane Surface(18) = {17};
//+
Line Loop(19) = {9, 12, 6, 7};
//+
Plane Surface(20) = {19};
//+
Physical Line("left") = {7, 8};
//+
Physical Line("right") = {3, 4};
//+
Physical Line("bottom") = {1, 2};
//+
Physical Line("top") = {5, 6};
//+
Physical Line("xmid") = {11, 12};
//+
Physical Line("ymid") = {9, 10};
//+
Physical Surface("beam") = {14, 16, 18, 20};
