
Include "streamers_2D.dat";

// Circle used to define the mechanical bar
//=========================================

If(Flag_Six_Code == 1)
	Point(1) = {x_0,  y_0,  0.0, lc_x_0};
	Point(2) = {x_0,  y_1,  0.0, lc_x_0};
	Point(3) = {x_1,  y_0,  0.0, lc_x_1};
	Point(4) = {x_1,  y_1,  0.0, lc_x_1};
	Point(5) = {x_2,  y_0,  0.0, lc_x_2};
	Point(6) = {x_2,  y_1,  0.0, lc_x_2};

	Line(1) = {1, 2};
	Line(2) = {1, 3};
	Line(3) = {2, 4};
	Line(4) = {3, 4};
	Line(5) = {3, 5};
	Line(6) = {4, 6};
	Line(7) = {5, 6};

	ll_internal = newll;
	Curve Loop(ll_internal) = {2, 4, -3, -1};
	surf_internal = news ; Plane Surface(news) = {ll_internal};

	ll_external = newll;
	Curve Loop(ll_external) = {5, 7, -6, -4};
	surf_external= news ; Plane Surface(news) = {ll_external};

	// Defining physical regions
	//==========================
	Physical Line(GAMMA_UP) = {3, 6};
	Physical Line(GAMMA_RIGHT) = {7};
	Physical Line(GAMMA_LEFT) = {1};
	Physical Line(GAMMA_BOTTOM) = {2, 5};
	Physical Surface(OMEGA_INTERNAL) = {surf_internal};
	Physical Surface(OMEGA_EXTERNAL) = {surf_external};
EndIf
If(Flag_Six_Code == 0)
	p = newp;
	divfactor = 5;
	flag_fine = 0;

	Point(p+1) = {x_0,  y_0,  0.0, lc_e_gas/divfactor};
	Point(p+2) = {x_2,  y_0,  0.0, lc_e_gas/divfactor};
	If(flag_fine == 0)
		Point(p+3) = {x_3,  y_0,  0.0, lc_e_external/(divfactor/4)};
	Else
		Point(p+3) = {x_3,  y_0,  0.0, lc_e_gas/(divfactor/1)};
	EndIf
	Point(p+4) = {x_4,  y_0,  0.0, lc_e_external/(divfactor/10)};

	Point(p+5) = {x_0,  y_1,  0.0, lc_e_gas/(divfactor)};
	Point(p+6) = {x_2,  y_1,  0.0, lc_e_gas/(divfactor)};

	Point(p+7) = {x_0,  y_2,  0.0, lc_e_gas/(divfactor)};
	Point(p+8) = {x_0,  y_3,  0.0, lc_e_gas/(divfactor)};
	Point(p+9) = {x_1,  y_3,  0.0, lc_e_gas/(divfactor)};

	If(flag_fine == 0)
		Point(p+10) = {x_1,  y_4,  0.0, lc_e_gas/(divfactor/1)};
		Point(p+11) = {x_3,  y_4,  0.0, lc_e_gas/(divfactor/1)};
		//Point(p+11) = {x_3,  y_4,  0.0, lc_e_external/(divfactor/2)};
	Else
		Point(p+10) = {x_1,  y_4,  0.0, lc_e_gas/(divfactor/1)};
		Point(p+11) = {x_3,  y_4,  0.0, lc_e_gas/(divfactor/1)};
	EndIf
	Point(p+12) = {x_4,  y_4,  0.0, lc_e_external};

	l = newl;
	Line(l+1) = {p+1, p+2};
	Line(l+2) = {p+2, p+3};
	Line(l+3) = {p+3, p+4};
	Line(l+4) = {p+1, p+5};
	Line(l+5) = {p+2, p+6};
	Line(l+6) = {p+5, p+6};
	Line(l+7) = {p+5, p+7};
	Line(l+8) = {p+3, p+11};
	Line(l+9) = {p+4, p+12};
	Circle(l+10) = {p+7, p+8, p+9};
	Line(l+11) = {p+9, p+10};
	Line(l+12) = {p+10, p+11};
	Line(l+13) = {p+11, p+12};

	Line Loop(newll) = {(l+1), (l+5), -(l+6), -(l+4)};
	surf_dielectrics = news ; Plane Surface(news) = {newll-1};

	Line Loop(newll) = {l+3, l+9, -(l+13), -(l+8)};
	surf_external= news ; Plane Surface(news) = {newll-1};

	Line Loop(newll) = {l+2, l+8, -(l+12), -(l+11), -(l+10), -(l+7), l+6, -(l+5)};
	surf_gas= news ; Plane Surface(news) = {newll-1};



	// Defining physical regions
	//==========================
	Physical Line(GAMMA_UP_1) = {l+10, l+11};
	Physical Line(GAMMA_UP_2) = {l+12, l+13};
	Physical Line(GAMMA_RIGHT) = {l+9};
	Physical Line(GAMMA_LEFT) = {l+4, l+7};
	Physical Line(GAMMA_BOTTOM) = {l+1, l+2, l+3};
	Physical Surface(OMEGA_DIELECTRICS) = {surf_dielectrics};
	Physical Surface(OMEGA_GAS) = {surf_gas};
	Physical Surface(OMEGA_EXTERNAL) = {surf_external};


	eps = 1.0e-6;
	pmesh_1[] += newp; Point(newp) = {0.5*x_2,  y_0 + eps,  0.0, lc_e_gas/(divfactor)};
	pmesh_1[] += newp; Point(newp) = {0.5*x_2,  0.5*(y_0+0.5*y_1) + eps,  0.0, lc_e_gas/(divfactor)};
	pmesh_1[] += newp; Point(newp) = {0.5*x_2,  0.5*y_1,  0.0, lc_e_gas/(divfactor)};
	pmesh_1[] += newp; Point(newp) = {0.5*x_2,  (y_1 - 0.25*y_1) - eps,  0.0, lc_e_gas/(divfactor)};
	pmesh_1[] += newp; Point(newp) = {0.5*x_2,  y_1 - eps,  0.0, lc_e_gas/(divfactor)};
	pmesh_2[] += newp; Point(newp) = {0.5*(x_1+x_2),  0.75*y_1 + 0.25*y_2,  0.0, lc_e_gas/(divfactor)};
	pmesh_2[] += newp; Point(newp) = {0.5*(x_1+x_2),  0.5*y_1+0.5*y_2,  0.0, lc_e_gas/(divfactor)};
	pmesh_2[] += newp; Point(newp) = {0.5*x_2,  y_2,  0.0, lc_e_gas/(divfactor/1)};

	Point{pmesh_1[]} In Surface {surf_dielectrics};
	Point{pmesh_2[]} In Surface {surf_gas};
EndIf
