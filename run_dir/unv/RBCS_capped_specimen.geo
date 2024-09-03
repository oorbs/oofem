//+
Point(1) = {40, -40, 0, 1.0};
//+
Point(2) = {60, -40, 0, 1.0};
//+
Point(3) = {100, 0, 0, 1.0};
//+
Point(4) = {0, 0, 0, 1.0};
//+
Point(5) = {0, 200, 0, 1.0};
//+
Point(6) = {100, 200, 0, 1.0};
//+
Point(7) = {60, 240, 0, 1.0};
//+
Point(8) = {40, 240, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {3, 6};
//+
Line(6) = {6, 5};
//+
Line(7) = {5, 4};
//+
Line(8) = {6, 7};
//+
Line(9) = {7, 8};
//+
Line(10) = {8, 5};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, -7, -6, -5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {6, -10, -9, -8};
//+
Plane Surface(3) = {3};
//+
Physical Curve("Top_ridge", 11) = {9};
//+
Physical Curve("Bot_ridge", 12) = {1};
//+
Physical Surface("Bot_cap", 13) = {1};
//+
Physical Surface("Specimen", 14) = {2};
//+
Physical Surface("Top_cap", 15) = {3};
