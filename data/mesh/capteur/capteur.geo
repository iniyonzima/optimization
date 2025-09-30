Include "capteur_data.geo";

Mesh.MshFileVersion = 2.2;

SetFactory("OpenCASCADE");

Geometry.CopyMeshingMethod = 1;

mm = 1.0e-3;

radius_source = 0.1 * mm; // radius of the primary inductor
radius_ext = 50 * mm; // radius of the disk that defines the entire domain

Mesh.CharacteristicLengthMax = radius_ext/10;

L_indu = 0.2 * mm; // lenght of the lateral inductors
W_indu = 0.4 * mm; // width of the lateral inductors

L_Fe = 5.0 * mm; // length of the ferromagnetic domain
W_Fe = 0.8 * mm; // width of the ferromagnetic domain

// Create the primary inductor
//============================
Disk_Int = news;
Disk(Disk_Int) = {0, -0, 0, radius_source};

// Create the secondary inductors
//===============================
Rect_R_L = news;
Rectangle(Rect_R_L) = {1.0 * mm, -0.2 * mm, 0, L_indu, W_indu, 0};
Rect_R_R = news;
Rectangle(Rect_R_R) = {1.4 * mm, -0.2 * mm, 0, L_indu, W_indu, 0};
Rect_L_R = news;
Rectangle(Rect_L_R) = {-1.2 * mm, -0.2 * mm, 0, L_indu, W_indu, 0};
Rect_L_L = news;
Rectangle(Rect_L_L) = {-1.6 * mm, -0.2 * mm, 0, L_indu, W_indu, 0};

// Create the ferromagnetic domains
//=================================
Rect_Fe_Bottom = news;
Rectangle(Rect_Fe_Bottom) = {-2.5 * mm, -1.2 * mm, 0, L_Fe, W_Fe, 0};
Rect_Fe_Top = news;
Rectangle(Rect_Fe_Top) = {-2.5 * mm, 0.4 * mm, 0, L_Fe, W_Fe, 0};

// Create the disk containing the surrounding air
//===============================================
Disk_Ext = newc;
Disk(Disk_Ext) = {0, -0, 0, radius_ext};

// Boolean operations to get the different subdomains
//=================================================== 
surf_all() = BooleanFragments{ Surface{Disk_Ext}; Delete; }{ Surface{Disk_Int, Rect_R_L, Rect_R_R, Rect_L_L, Rect_L_R, Rect_Fe_Bottom, Rect_Fe_Top}; Delete; };

bnd_sphere_ext() = CombinedBoundary{ Surface{surf_all()}; };

bnd_primary_indu() = CombinedBoundary{ Surface{surf_all(0)}; };
MeshSize{ PointsOf{ Line{bnd_primary_indu()}; } } = radius_source/5;


pdfTeX warning (dest): name{Hfootnote.1} has been referenced but does not exist
, replaced by a fixed one

</usr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-Bold.pfb></u
sr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-BoldItalic.pfb>
</usr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-Italic.pfb><
/usr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-Regular.pfb><


pdfTeX warning (dest): name{Hfootnote.1} has been referenced but does not exist
, replaced by a fixed one

</usr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-Bold.pfb></u
sr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-BoldItalic.pfb>
</usr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-Italic.pfb><
/usr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-Regular.pfb><


pdfTeX warning (dest): name{Hfootnote.1} has been referenced but does not exist
, replaced by a fixed one

</usr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-Bold.pfb></u
sr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-BoldItalic.pfb>
</usr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-Italic.pfb><
/usr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-Regular.pfb><


pdfTeX warning (dest): name{Hfootnote.1} has been referenced but does not exist
, replaced by a fixed one

</usr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-Bold.pfb></u
sr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-BoldItalic.pfb>
</usr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-Italic.pfb><
/usr/share/texlive/texmf-dist/fonts/type1/public/stix/STIXGeneral-Regular.pfb><

bnd_secondary_indu() = CombinedBoundary{ Surface{surf_all(1), surf_all(2), surf_all(3), surf_all(4)}; };
MeshSize{ PointsOf{ Line{bnd_secondary_indu()}; } } = L_indu/(4 * 5);

bnd_fe() = CombinedBoundary{ Surface{surf_all(5), surf_all(6)}; };
MeshSize{ PointsOf{ Line{bnd_fe()}; } } = W_Fe/32;

MeshSize{ PointsOf{ Line{bnd_sphere_ext()}; } } = radius_ext/2;


// Definition of physical groups
//==============================
OMEGA_INDU_P = 2000;
OMEGA_INDU_S_L_L = 2001;
OMEGA_INDU_S_L_R = 2002;
OMEGA_INDU_S_R_L = 2003;
OMEGA_INDU_S_R_R = 2004;
OMEGA_MAG_TOP = 2005;
OMEGA_MAG_BOTTOM = 2006;
OMEGA_AIR = 2007;
GAMMA_INF = 1000;

Physical Surface('Primary inductor', OMEGA_INDU_P) = {surf_all(0)};
Physical Surface('Secondary inductor R L', OMEGA_INDU_S_R_L) = {surf_all(1)};
Physical Surface('Secondary inductor R R', OMEGA_INDU_S_R_R) = {surf_all(2)};
Physical Surface('Secondary inductor L L', OMEGA_INDU_S_L_L) = {surf_all(3)};
Physical Surface('Secondary inductor L R', OMEGA_INDU_S_L_R) = {surf_all(4)};
Physical Surface(OMEGA_MAG_BOTTOM) = {surf_all(5)};
Physical Surface(OMEGA_MAG_TOP) = {surf_all(6)};
Physical Surface(OMEGA_AIR) = {surf_all(7)};

Physical Line(GAMMA_INF) = {bnd_sphere_ext()};
