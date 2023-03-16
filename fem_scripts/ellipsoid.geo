SetFactory("OpenCASCADE");

DefineConstant[
  Sx = {4, Min 2, Max 6, Step 0.05, Name "Geometry/Box/Sx"},
  Sy = {4, Min 2, Max 6, Step 0.05, Name "Geometry/Box/Sy"},
  Sz = {4, Min 2, Max 6, Step 0.05, Name "Geometry/Box/Sz"},
  lcoarse = {0.5, Min 0.01, Max 1, Step 0.05, Name "Geometry/Parameters/lcoarse"},
  per2pi = {20, Min 10, Max 150, Step 5, Name "Geometry/Parameters/per2pi"},
  Rx = {1, Min 0.25, Max 1.5, Step 0.05, Name "Geometry/Inclusion/Rx"},
  Ry = {1, Min 0.25, Max 1.5, Step 0.05, Name "Geometry/Inclusion/Ry"},
  Rz = {1, Min 0.25, Max 1.5, Step 0.05, Name "Geometry/Inclusion/Rz"}
  tinterface = {0.5, Min 0.05, Max 2, Step 0.05, Name "Geometry/Inclusion/tinterface"}
];

Mesh.MshFileVersion = 2.2;

Mesh.Algorithm3D = 4;  // Frontal algorithm required for NetGen optimization.
// Mesh size is completely determined by background field.
// Mesh.CharacteristicLengthExtendFromBoundary = 1;
// Mesh.CharacteristicLengthFromPoints = 1;
Mesh.MinimumElementsPerTwoPi = per2pi;
Mesh.CharacteristicLengthFromCurvature = 1;

Point(1) = {0,0,0};
Box(1) = {-Sx/2,-Sy/2,-Sz/2, Sx/2,Sy/2,Sz/2};

// Boolean operations with OpenCASCADE always create new entities. Adding
// `Delete' in the arguments allows to automatically delete the original
// entities.

// Inclusion
Sphere(2) = {0, 0, 0, 1};
Dilate {{0, 0, 0}, {Rx, Ry, Rz}} { Volume{2}; }

// Interface
Sphere(3) = {0, 0, 0, 1};
Dilate {{0, 0, 0}, {Rx + tinterface, Ry + tinterface, Rz + tinterface}} { Volume{3}; }

// boxes used to cut out 1/8 of the cell
Box(4) = {0,-Sy,-Sz, 2*Sx,2*Sy,2*Sz};  // cut out the positive X
Box(5) = {-Sx,0,-Sz, 2*Sx,2*Sy,2*Sz};  // cut out the positive Y
Box(6) = {-Sx,-Sy,0, 2*Sx,2*Sy,2*Sz};  // cut out the positive Z
Box(7) = {-Sx, -Sy,-Sz, Sx/2,2*Sy,2*Sz};  // cut out the positive X
Box(8) = {-Sx,-Sy,-Sz, 2*Sx,Sy/2,2*Sz};  // cut out the positive Y
Box(9) = {-Sx,-Sy,-Sz, 2*Sx,2*Sy,Sz/2};  // cut out the positive Z


// Cut the sphere and combine the result with the box
v() = BooleanDifference{ Volume{2:3}; Delete; }{ Volume{4:9}; Delete; };
v() = BooleanFragments{ Volume{1:2}; Delete; }{ Volume{3}; Delete; };

// Assign a mesh size to all the points of all the volumes:
Characteristic Length{ PointsOf{ Volume{:}; } } = lcoarse;

Physical Volume("Matrix") = {4};
Physical Volume("Interface") = {3};
Physical Volume("Inclusion") = {2};

eps = 1e-2;
Sxmin() = Surface In BoundingBox{-Sx/2-eps, -Sy/2-eps, -Sz/2-eps, -Sx/2+eps, Sy/2+eps, Sz/2+eps};
Sxmax() = Surface In BoundingBox{0-eps, -Sy/2-eps, -Sz/2-eps, 0+eps, Sy/2+eps, Sz/2+eps};
Symin() = Surface In BoundingBox{-Sx/2-eps, -Sy/2-eps, -Sz/2-eps, Sx/2+eps, -Sy/2+eps, Sz/2+eps};
Symax() = Surface In BoundingBox{-Sx/2-eps, 0-eps, -Sz/2-eps, Sx/2+eps, 0+eps, Sz/2+eps};
Szmin() = Surface In BoundingBox{-Sx/2-eps, -Sy/2-eps, -Sz/2-eps, Sx/2+eps, Sy/2+eps, -Sz/2+eps};
Szmax() = Surface In BoundingBox{-Sx/2-eps, -Sy/2-eps, 0-eps, Sx/2+eps, Sy/2+eps, 0+eps};

//Printf("Bottom surface: %g %g; Top surface: %g %g %g", Szmin(0), Szmin(1), Szmax(0), Szmax(1), Szmax(2));
//Printf("Sides surface: %g %g %g %g %g %g %g %g %g %g", Sxmin(0), Sxmin(1), Sxmax(0), Sxmax(1), Sxmax(2), Symin(0), Symin(1), Symax(0), Symax(1), Symax(2));
sides() = {};
sides() += Sxmin();
sides() += Sxmax();
sides() += Symin();
sides() += Symax();
Physical Surface("Sides") = sides();
Physical Surface("Top") = Szmax();
Physical Surface("Bottom") = Szmin();


// Print some info about the model being used
Printf("%%%% START PARAMS");
Printf("inc_type = 'ellipsoid';");
Printf("order = 1;");
Printf("Sx = %g;", Sx);
Printf("Sy = %g;", Sy);
Printf("Sz = %g;", Sz);
Printf("lcoarse = %g;", lcoarse);
Printf("per2pi = %g;", per2pi);
Printf("Rx = %g;", Rx);
Printf("Ry = %g;", Ry);
Printf("Rz = %g;", Rz);
Printf("ninc = 1;");
Printf("ncells = 2;");
Printf("vl = %g;", 4/3*Pi*Rx*Ry*Rz / (Sx*Sy*Sz));
Printf("%%%% END PARAMS");
