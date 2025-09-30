
Include "streamers_2D_six_codes.dat";

Mesh.MshFileVersion = 2.2;

Num_Transfinite_Y_Bottom = 4000/1.3; // 2500;
Num_Transfinite_Y_Top = 1000/1.3; // 2500;
Num_Transfinite_X = 5000; // 2500;


// Circle used to define the mechanical bar
//=========================================

Point(1) = {x_0,  y_0,  0.0, lc_x_0};
Point(2) = {x_0,  y_1,  0.0, lc_x_0};
Point(3) = {x_1,  y_0,  0.0, lc_x_1};
Point(4) = {x_1,  y_1,  0.0, lc_x_1};
Point(5) = {x_2,  y_0,  0.0, lc_x_2};
Point(6) = {x_2,  y_1,  0.0, lc_x_2};

If(Flag_Refined_Geometry)
	Point(7) = {x_0,  y_2,  0.0, lc_x_0};
	Point(8) = {x_1,  y_2,  0.0, lc_x_1};
	Point(9) = {x_2,  y_2,  0.0, lc_x_2};
EndIf

Line(1) = {1, 2};
Line(2) = {1, 3};
Line(3) = {2, 4};
Line(4) = {3, 4};
Line(5) = {3, 5};
Line(6) = {4, 6};
Line(7) = {5, 6};
If(Flag_Refined_Geometry)
	Line(8) = {2, 7};
	Line(9) = {7, 8};
	Line(10) = {4, 8};
	Line(11) = {8, 9};
	Line(12) = {6, 9};
EndIf



If(Flag_Refined_Geometry)
	ll_internal = newll;
	Curve Loop(ll_internal) = {2, 4, -3, -1};
	surf_internal = news ; Plane Surface(news) = {ll_internal};
	
	If(Flag_Transfinite)
		Transfinite Curve {1, 4} = Num_Transfinite_Y_Bottom Using Progression 0.9991956047;
		Transfinite Curve {2, 3} = Floor[ (Num_Transfinite_X/(2 * Ratio_Reg1_Reg2) ) ] Using Progression 1.009;
		Transfinite Surface {news-1};
		Recombine Surface {news-1};
	EndIf

	ll_internal_other = newll;
	Curve Loop(ll_internal_other) = {3, 10, -9, -8};
	surf_internal_other = news ; Plane Surface(news) = {ll_internal_other};		
	If(Flag_Transfinite)
		Transfinite Curve {8, 10} = Num_Transfinite_Y_Top Using Progression 1.0032240619;
		Transfinite Curve {9} = (Num_Transfinite_X/(2 * Ratio_Reg1_Reg2)) Using Progression 1.009;
		Transfinite Surface {news-1};
		Recombine Surface {news-1};
	EndIf

	ll_external = newll;
	Curve Loop(ll_external) = {5, 7, -6, -4};
	surf_external= news ; Plane Surface(news) = {ll_external};

	ll_external_other = newll;
	Curve Loop(ll_external_other) = {6, 12, -11, -10};
	surf_external_other = news ; Plane Surface(news) = {ll_external_other};
EndIf

If(Flag_Refined_Geometry == 0)
	ll_internal = newll;
	Curve Loop(ll_internal) = {2, 4, -3, -1};
	surf_internal = news ; Plane Surface(news) = {ll_internal};
	
	ll_external = newll;
	Curve Loop(ll_external) = {5, 7, -6, -4};
	surf_external= news ; Plane Surface(news) = {ll_external};
	
	//Transfinite Curve {1, 4, 7} = (60*5) Using Progression 1.0;
	//Transfinite Curve {5, 6} = (30) Using Progression 1.0;
	//Transfinite Curve {2, 3} = (30) Using Progression 1.0;
	Transfinite Curve {1, 4} = (num_Y) Using Progression 1.0;
	Transfinite Curve {5, 6} = (num_X_2) Using Progression 1.0;
	Transfinite Curve {2, 3} = (num_X_1) Using Progression 1.0;
	//Transfinite Surface {surf_internal, surf_external};
	If(Flag_Transfinite)
		Recombine Surface {surf_internal, surf_external};
	EndIf
EndIf

///*

// Defining physical regions
//==========================
If(Flag_Refined_Geometry == 0)
	Physical Surface(OMEGA_INTERNAL) = {surf_internal};
	Physical Surface(OMEGA_EXTERNAL) = {surf_external};
	Physical Line(GAMMA_UP) = {3, 6};
	Physical Line(GAMMA_RIGHT) = {7};
	Physical Line(GAMMA_LEFT) = {1};
	Physical Line(GAMMA_BOTTOM) = {2, 5};
EndIf
If(Flag_Refined_Geometry)
	Physical Surface(OMEGA_INTERNAL) = {surf_internal, surf_internal_other};
	Physical Surface(OMEGA_EXTERNAL) = {surf_external, surf_external_other};
	Physical Line(GAMMA_UP) = {9, 11};
	Physical Line(GAMMA_RIGHT) = {7, 12};
	Physical Line(GAMMA_LEFT) = {1, 8};
	Physical Line(GAMMA_BOTTOM) = {2, 5};
EndIf

//*/