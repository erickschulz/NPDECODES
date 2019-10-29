//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0.0, 0.0, 0, 1, 1, 0};
//+
Physical Curve("neu") = {3, 2};
//+
Physical Curve("neu") -= {2};
//+
Physical Surface("area") = {1};
//+
Physical Curve("dir") = {4, 1, 2};
