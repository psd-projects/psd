// Plane-strain strip-footing benchmark of Cermak, Sysala, and Valdman.
// The domain is the right half of a symmetric footing problem.
//
// Boundary labels shared by DOLFINx and PSD:
//   1: bottom, uy = 0
//   2: right,  ux = 0
//   3: footing, uy = prescribed settlement
//   4: traction-free top
//   5: left symmetry, ux = 0
//   6: material domain

SetFactory("Built-in");

domain_size = 10.0;
footing_width = 1.0;
cells_per_side = 10;

Point(1) = {0.0,           0.0,         0.0};
Point(2) = {domain_size,   0.0,         0.0};
Point(3) = {domain_size,   domain_size, 0.0};
Point(4) = {footing_width, domain_size, 0.0};
Point(5) = {0.0,           domain_size, 0.0};

Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(4) = {3, 4}; // free top
Line(3) = {4, 5}; // footing
Line(5) = {5, 1}; // left

Curve Loop(1) = {1, 2, 4, 3, 5};
Plane Surface(6) = {1};

// Ten intervals on each full-height/width boundary. The split top uses
// one interval below the footing and nine along the free surface.
Transfinite Curve {1} = cells_per_side + 1;
Transfinite Curve {2, 5} = cells_per_side + 1;
Transfinite Curve {3} = cells_per_side/10 + 1;
Transfinite Curve {4} = 9*cells_per_side/10 + 1;
Transfinite Surface {6} = {1, 2, 3, 5};

Physical Curve("bottom", 1) = {1};
Physical Curve("right", 2) = {2};
Physical Curve("footing", 3) = {3};
Physical Curve("top_free", 4) = {4};
Physical Curve("left", 5) = {5};
Physical Surface("domain", 6) = {6};

Mesh.Algorithm = 8;
Mesh.RecombineAll = 0;
Mesh.MshFileVersion = 2.2;
