// French Alps mesh with a conformal basin, using only the Gmsh .geo language.
//
// The original 200 x 200 triangulated terrain has been sampled onto a 21 x 21
// grid. The samples below are baked into this file so that every entity is an
// OpenCASCADE CAD entity and BooleanFragments can create two true volumes.
//
// Run with:
//   gmsh french_alps_two_volumes_boolean.geo -

SetFactory("OpenCASCADE");
// Use tight triangulation-based bounds for trimmed Boolean faces.
Geometry.OCCBoundsUseStl = 1;

nx = 21;
ny = 21;
xMin = -0.198781;
xMax =  3.29848;
yMin = -3.3514;
yMax =  0.145877;
zBottom = -1.0;

dx = (xMax - xMin) / (nx - 1);
dy = (yMax - yMin) / (ny - 1);

// Row-major terrain elevations: index = i + nx * j.
zSamples[] = {
  0.200829, 0.208141436, 0.245470278, 0.28879704, 0.319999095, 0.309437559, 0.297237983, 0.315652219, 0.348404651, 0.344692703, 0.311490142, 0.244679352, 0.186090101, 0.156391843, 0.139492472, 0.1565527, 0.202542488, 0.222485876, 0.227144272, 0.232327859, 0.22776,
  0.191855705, 0.180225607, 0.233645226, 0.255296611, 0.295185667, 0.327476839, 0.322687425, 0.32678001, 0.369740342, 0.394045363, 0.363204895, 0.262344979, 0.207858717, 0.168186178, 0.159571575, 0.141092527, 0.195776266, 0.232553744, 0.233685127, 0.224473145, 0.22344866,
  0.179707428, 0.178727655, 0.162764601, 0.206279625, 0.264923921, 0.325737224, 0.369307709, 0.356178462, 0.448739122, 0.427434837, 0.326227584, 0.293765854, 0.293876267, 0.242303436, 0.200928012, 0.165593877, 0.159927576, 0.172048006, 0.214259518, 0.223903854, 0.22202272,
  0.163589199, 0.155685036, 0.16725828, 0.246200269, 0.294425295, 0.314093754, 0.33136621, 0.426030174, 0.400883035, 0.431612897, 0.353835339, 0.360292708, 0.323535379, 0.308585588, 0.269451773, 0.237243615, 0.200851057, 0.177066261, 0.163085515, 0.212044946, 0.217075134,
  0.154273164, 0.165883126, 0.213309731, 0.199128652, 0.223386332, 0.286514585, 0.34972558, 0.367382135, 0.364620855, 0.393143242, 0.389654187, 0.335331513, 0.319499979, 0.331031617, 0.316466139, 0.32139846, 0.272200414, 0.236826048, 0.183020195, 0.171486276, 0.208272397,
  0.155813386, 0.13088809, 0.136430567, 0.154594981, 0.209038127, 0.263634272, 0.270737437, 0.288518685, 0.309060931, 0.342921436, 0.363282424, 0.328940227, 0.29355255, 0.318393518, 0.345616512, 0.391022718, 0.377625987, 0.334123171, 0.242144021, 0.203168656, 0.223456457,
  0.173942661, 0.173908935, 0.164897714, 0.166770967, 0.169447265, 0.159771131, 0.174724556, 0.236031845, 0.258536183, 0.276724455, 0.35564251, 0.32682813, 0.276358808, 0.255415661, 0.282609766, 0.326181717, 0.294282166, 0.338016794, 0.287236242, 0.267932393, 0.252009626,
  0.191092193, 0.183738321, 0.152579006, 0.121787331, 0.109395637, 0.114602793, 0.134742102, 0.14382941, 0.191043712, 0.239337671, 0.273257365, 0.34857825, 0.296627657, 0.225255712, 0.272927629, 0.286921173, 0.253513265, 0.300234435, 0.314967057, 0.279244075, 0.264339621,
  0.184347378, 0.183564184, 0.124997493, 0.101753467, 0.116064891, 0.0985669532, 0.101264796, 0.112922766, 0.13747368, 0.189500121, 0.237913611, 0.281638581, 0.298621919, 0.216198769, 0.219395049, 0.224648644, 0.270856616, 0.330456807, 0.298125904, 0.270114452, 0.278325483,
  0.16224878, 0.152117326, 0.101711641, 0.103021613, 0.138889582, 0.160209343, 0.136717438, 0.101605656, 0.109507814, 0.155485204, 0.204623215, 0.244247328, 0.26829368, 0.208631921, 0.23923539, 0.257986032, 0.271891676, 0.29313368, 0.352568213, 0.305286462, 0.30863538,
  0.134223381, 0.120957858, 0.0819612506, 0.125106302, 0.182840603, 0.22373154, 0.199969103, 0.145084997, 0.107438921, 0.103954052, 0.166271263, 0.211308494, 0.200465062, 0.194531911, 0.256668858, 0.295869953, 0.29500617, 0.300265246, 0.352745479, 0.346866845, 0.327358636,
  0.110284528, 0.0835034167, 0.0806881224, 0.107019957, 0.151620601, 0.197314522, 0.212950168, 0.198583321, 0.150741264, 0.104069435, 0.114244117, 0.182567225, 0.176313539, 0.222476589, 0.320217387, 0.368853061, 0.381347586, 0.35043277, 0.2963501, 0.296209059, 0.330639233,
  0.11605518, 0.103946232, 0.0978429742, 0.139032634, 0.169834255, 0.136743434, 0.182929353, 0.237606512, 0.185590325, 0.12622415, 0.106421122, 0.113263429, 0.149869023, 0.219834554, 0.273136455, 0.348897668, 0.291037054, 0.270689588, 0.279596939, 0.326918913, 0.33695238,
  0.144012937, 0.140637989, 0.143271258, 0.160251471, 0.210793597, 0.193444279, 0.160318259, 0.220182697, 0.220396847, 0.17187962, 0.12869688, 0.122181457, 0.19532761, 0.239986393, 0.276956045, 0.290257979, 0.255048796, 0.290762032, 0.316839623, 0.339835904, 0.337342317,
  0.172791529, 0.189382234, 0.182959603, 0.171313537, 0.195377989, 0.201781376, 0.172687332, 0.227744063, 0.258567452, 0.215217778, 0.182545336, 0.124647821, 0.15858726, 0.207325316, 0.230967601, 0.234061241, 0.255732119, 0.313717339, 0.359340575, 0.318790863, 0.319102067,
  0.189583831, 0.201060171, 0.215945108, 0.262100092, 0.200405913, 0.19797297, 0.179288187, 0.225135532, 0.266620667, 0.252639926, 0.219518927, 0.170106374, 0.120926482, 0.14573664, 0.176743372, 0.21467865, 0.286201801, 0.323572811, 0.328217687, 0.303049811, 0.307019719,
  0.196103848, 0.19289786, 0.221849671, 0.259755752, 0.220390482, 0.231482727, 0.201244172, 0.206644345, 0.230538752, 0.275720186, 0.24459092, 0.214674003, 0.153586071, 0.140307857, 0.211149358, 0.265826614, 0.283663465, 0.299330557, 0.32369688, 0.326209424, 0.30371398,
  0.201645613, 0.193340465, 0.223974789, 0.210783091, 0.209329839, 0.227849834, 0.207446019, 0.218863561, 0.241881274, 0.228478133, 0.266284567, 0.229434927, 0.157742425, 0.136716634, 0.175434406, 0.20990206, 0.264035589, 0.297590885, 0.336692387, 0.305881878, 0.297304666,
  0.204228092, 0.195381823, 0.1887723, 0.181431582, 0.190969739, 0.198708243, 0.212846353, 0.216968791, 0.242182809, 0.187185847, 0.209526624, 0.223168346, 0.148768607, 0.178878636, 0.160593225, 0.224577139, 0.28817676, 0.283031435, 0.301696397, 0.280943846, 0.276441934,
  0.19729516, 0.201730592, 0.125839327, 0.172222804, 0.15733541, 0.143753219, 0.198755762, 0.273970363, 0.266249502, 0.254194954, 0.197972435, 0.155812303, 0.134998469, 0.201074677, 0.184322526, 0.212838857, 0.256401452, 0.225275946, 0.218700661, 0.228275464, 0.242977013,
  0.18336, 0.16996653, 0.130872876, 0.114135263, 0.118213349, 0.139319204, 0.185835423, 0.245953004, 0.267040368, 0.252693782, 0.211420179, 0.161845998, 0.148719863, 0.177464067, 0.198592664, 0.214662965, 0.237049633, 0.223408915, 0.196696482, 0.204052591, 0.222291
};

// Build a shared-edge triangulated CAD surface.
For j In {0:ny - 1}
  For i In {0:nx - 1}
    k = i + nx * j;
    p = newp;
    Point(p) = {xMin + i * dx, yMin + j * dy, zSamples[k]};
    topPoints[k] = p;
  EndFor
EndFor

For j In {0:ny - 1}
  For i In {0:nx - 2}
    k = i + (nx - 1) * j;
    l = newl;
    Line(l) = {topPoints[i + nx * j], topPoints[i + 1 + nx * j]};
    xEdges[k] = l;
  EndFor
EndFor

For j In {0:ny - 2}
  For i In {0:nx - 1}
    k = i + nx * j;
    l = newl;
    Line(l) = {topPoints[i + nx * j], topPoints[i + nx * (j + 1)]};
    yEdges[k] = l;
  EndFor
EndFor

For j In {0:ny - 2}
  For i In {0:nx - 2}
    k = i + (nx - 1) * j;
    l = newl;
    Line(l) = {topPoints[i + nx * j], topPoints[i + 1 + nx * (j + 1)]};
    diagonalEdges[k] = l;
  EndFor
EndFor

For j In {0:ny - 2}
  For i In {0:nx - 2}
    cell = i + (nx - 1) * j;

    loopA = newll;
    Curve Loop(loopA) = {
      xEdges[i + (nx - 1) * j],
      yEdges[i + 1 + nx * j],
      -diagonalEdges[cell]
    };
    sA = news;
    Plane Surface(sA) = {loopA};
    topSurfaces[2 * cell] = sA;

    loopB = newll;
    Curve Loop(loopB) = {
      diagonalEdges[cell],
      -xEdges[i + (nx - 1) * (j + 1)],
      -yEdges[i + nx * j]
    };
    sB = news;
    Plane Surface(sB) = {loopB};
    topSurfaces[2 * cell + 1] = sB;
  EndFor
EndFor

// Bottom corners, vertical corner curves and the rectangular bottom.
baseSW = newp; Point(baseSW) = {xMin, yMin, zBottom};
baseSE = newp; Point(baseSE) = {xMax, yMin, zBottom};
baseNE = newp; Point(baseNE) = {xMax, yMax, zBottom};
baseNW = newp; Point(baseNW) = {xMin, yMax, zBottom};

bSouth = newl; Line(bSouth) = {baseSW, baseSE};
bEast  = newl; Line(bEast)  = {baseSE, baseNE};
bNorth = newl; Line(bNorth) = {baseNE, baseNW};
bWest  = newl; Line(bWest)  = {baseNW, baseSW};

vSW = newl; Line(vSW) = {baseSW, topPoints[0]};
vSE = newl; Line(vSE) = {baseSE, topPoints[nx - 1]};
vNE = newl; Line(vNE) = {baseNE, topPoints[nx * ny - 1]};
vNW = newl; Line(vNW) = {baseNW, topPoints[nx * (ny - 1)]};

bottomLoop = newll;
Curve Loop(bottomLoop) = {bSouth, bEast, bNorth, bWest};
bottomSurface = news;
Plane Surface(bottomSurface) = {bottomLoop};

southCurves[] = {bSouth, vSE};
For i In {nx - 2:0:-1}
  southCurves[] += {-xEdges[i]};
EndFor
southCurves[] += {-vSW};
southLoop = newll; Curve Loop(southLoop) = {southCurves[]};
southSurface = news; Plane Surface(southSurface) = {southLoop};

eastCurves[] = {bEast, vNE};
For j In {ny - 2:0:-1}
  eastCurves[] += {-yEdges[nx - 1 + nx * j]};
EndFor
eastCurves[] += {-vSE};
eastLoop = newll; Curve Loop(eastLoop) = {eastCurves[]};
eastSurface = news; Plane Surface(eastSurface) = {eastLoop};

northCurves[] = {bNorth, vNW};
For i In {0:nx - 2}
  northCurves[] += {xEdges[i + (nx - 1) * (ny - 1)]};
EndFor
northCurves[] += {-vNE};
northLoop = newll; Curve Loop(northLoop) = {northCurves[]};
northSurface = news; Plane Surface(northSurface) = {northLoop};

westCurves[] = {bWest, vSW};
For j In {0:ny - 2}
  westCurves[] += {yEdges[nx * j]};
EndFor
westCurves[] += {-vNW};
westLoop = newll; Curve Loop(westLoop) = {westCurves[]};
westSurface = news; Plane Surface(westSurface) = {westLoop};

alpsShell = newsl;
Surface Loop(alpsShell) = {
  topSurfaces[], bottomSurface,
  southSurface, eastSurface, northSurface, westSurface
} Using Sewing;
alpsVolume = newv;
Volume(alpsVolume) = {alpsShell};

// Basin geometry from the Python version. The basin tool is the union of a
// lower-half ellipsoid and an elliptical cylinder extending above the terrain.
basinCenterX = 1.55;
basinCenterY = -1.60;
basinRadiusX = 1.15;
basinRadiusY = 0.95;
basinRimZ = 0.06;
basinDepth = 0.48;
toolTop = 0.65;

lowerEllipsoid = newv;
Sphere(lowerEllipsoid) = {
  basinCenterX, basinCenterY, basinRimZ, 1,
  -Pi / 2, 0, 2 * Pi
};
Dilate {
  {basinCenterX, basinCenterY, basinRimZ},
  {basinRadiusX, basinRadiusY, basinDepth}
} { Volume{lowerEllipsoid}; }

upperCylinder = newv;
Cylinder(upperCylinder) = {
  basinCenterX, basinCenterY, basinRimZ,
  0, 0, toolTop - basinRimZ, 1, 2 * Pi
};
Dilate {
  {basinCenterX, basinCenterY, basinRimZ},
  {basinRadiusX, basinRadiusY, 1}
} { Volume{upperCylinder}; }

basinTool[] = BooleanUnion {
  Volume{lowerEllipsoid}; Delete;
}{
  Volume{upperCylinder}; Delete;
};

// Fragmenting produces the background, the part of the basin tool inside the
// Alps, and the unused portion of the tool above the terrain.
fragments[] = BooleanFragments {
  Volume{alpsVolume}; Delete;
}{
  Volume{basinTool[]}; Delete;
};
Printf("BooleanFragments returned %g volumes", #fragments[]);

// Discard the portion of the tool above the terrain. It is the only fragment
// whose upper bounding box reaches toolTop.
domainVolumes[] = {};
For i In {0:#fragments[] - 1}
  bbox[] = BoundingBox Volume{fragments[i]};
  If(bbox[5] > 0.5)
    Recursive Delete { Volume{fragments[i]}; }
  Else
    domainVolumes[] += {fragments[i]};
  EndIf
EndFor

If(#domainVolumes[] != 2)
  Error("Boolean construction did not produce exactly two Alps volumes");
EndIf

// The basin is wholly inside its elliptical X/Y bounding box; the background
// reaches the rectangular Alps boundary.
basinVolumes[] = {};
backgroundVolumes[] = {};
For i In {0:#domainVolumes[] - 1}
  bbox[] = BoundingBox Volume{domainVolumes[i]};
  If(bbox[0] > xMin + 0.1 && bbox[3] < xMax - 0.1 &&
     bbox[1] > yMin + 0.1 && bbox[4] < yMax - 0.1)
    basinVolumes[] += {domainVolumes[i]};
  Else
    backgroundVolumes[] += {domainVolumes[i]};
  EndIf
EndFor

If(#basinVolumes[] != 1 || #backgroundVolumes[] != 1)
  Error("Could not identify the basin and background Boolean fragments");
EndIf

Physical Volume("volume_1_background", 11) = {backgroundVolumes[]};
Physical Volume("volume_2_basin", 22) = {basinVolumes[]};

// Recover external boundary groups from the finished Boolean model.
eps = 1e-6;
outerSurfaces[] = Abs(CombinedBoundary { Volume{domainVolumes[]}; });
bottomSurfaces[] = Surface In BoundingBox {
  xMin - eps, yMin - eps, zBottom - eps,
  xMax + eps, yMax + eps, zBottom + eps
};
southSurfaces[] = Surface In BoundingBox {
  xMin - eps, yMin - eps, zBottom - eps,
  xMax + eps, yMin + eps, toolTop + eps
};
eastSurfaces[] = Surface In BoundingBox {
  xMax - eps, yMin - eps, zBottom - eps,
  xMax + eps, yMax + eps, toolTop + eps
};
northSurfaces[] = Surface In BoundingBox {
  xMin - eps, yMax - eps, zBottom - eps,
  xMax + eps, yMax + eps, toolTop + eps
};
westSurfaces[] = Surface In BoundingBox {
  xMin - eps, yMin - eps, zBottom - eps,
  xMin + eps, yMax + eps, toolTop + eps
};

paraxialSurfaces[] = {
  bottomSurfaces[], southSurfaces[], eastSurfaces[],
  northSurfaces[], westSurfaces[]
};
topBoundarySurfaces[] = outerSurfaces[];
topBoundarySurfaces[] -= paraxialSurfaces[];

Physical Surface("paraxial", 1) = {paraxialSurfaces[]};
Physical Surface("top", 2) = {topBoundarySurfaces[]};

basinBoundary[] = Abs(Boundary { Volume{basinVolumes[]}; });
basinInterface[] = basinBoundary[];
basinInterface[] -= outerSurfaces[];
Physical Surface("basin_interface", 3) = {basinInterface[]};

// Double-couple source points, all located in the background volume.
dcLength = 0.01;
sourceX = 0.0;
sourceY = -2.8;
sourceZ = -0.8;

dcCenter = newp; Point(dcCenter) = {sourceX, sourceY, sourceZ, dcLength};
dcNorth = newp; Point(dcNorth) = {sourceX, sourceY, sourceZ + dcLength, dcLength};
dcSouth = newp; Point(dcSouth) = {sourceX, sourceY, sourceZ - dcLength, dcLength};
dcEast = newp; Point(dcEast) = {sourceX - dcLength, sourceY, sourceZ, dcLength};
dcWest = newp; Point(dcWest) = {sourceX + dcLength, sourceY, sourceZ, dcLength};

Physical Point("center", 2) = {dcCenter};
Physical Point("dcNorth", 3) = {dcNorth};
Physical Point("dcSouth", 4) = {dcSouth};
Physical Point("dcEast", 5) = {dcEast};
Physical Point("dcWest", 6) = {dcWest};

Point {dcCenter, dcNorth, dcSouth, dcEast, dcWest}
  In Volume {backgroundVolumes[0]};

Mesh.MeshSizeMin = 0.04;
Mesh.MeshSizeMax = 0.16;
General.NumThreads = 8;
Mesh.MshFileVersion = 2.2;
Mesh.Binary = 1;
Mesh.Algorithm = 5;
Mesh.Algorithm3D = 10;

Mesh 3;
Save "french_alps_two_volumes.msh";

Printf("Boolean volume tags: background=%g, basin=%g",
       backgroundVolumes[0], basinVolumes[0]);
