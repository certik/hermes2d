#include "forms.cpp"
/*** Boundary types and conditions ***/

// Poisson takes Dirichlet and Neumann boundaries
int phi_bc_types(int marker) {
	return (marker == SIDE_MARKER || marker == TOP_MARKER) 
		? BC_NATURAL : BC_ESSENTIAL;
}

// Nernst-Planck takes Neumann boundaries
int C_bc_types(int marker) {
	return BC_NATURAL;
}

// Diricleht Boundary conditions for Poisson equation.
scalar phi_bc_values(int marker, double x, double y) {
	return 0.0;
}

template<class Real, class Scalar>
Scalar linear_form_surf_top(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
	return -E_FIELD * int_v<Real, Scalar>(n, wt, v);
}

void solveNonadaptive(Mesh &mesh, NonlinSystem &nls,
		Solution &Cp, Solution &Ci, Solution &phip, Solution &phii) {
	begin_time();

	//VectorView vview("electric field [V/m]", 0, 0, 600, 600);
	ScalarView Cview("Concentration [mol/m3]", 0, 0, 800, 800);
	ScalarView phiview("Voltage [V]", 650, 0, 600, 600);

	#ifdef CONT_OUTPUT
	phiview.show(&phii);
	Cview.show(&Ci);
	Cview.wait_for_keypress();
	#endif
	Solution Csln, phisln;
	for (int n = 1; n <= NSTEP; n++) {

		#ifdef VERBAL
		info("\n---- Time step %d ----", n);
		#endif

		int it = 1;
		double res_l2_norm;

		do {

			#ifdef VERBAL
			info("\n -------- Time step %d, Newton iter %d --------\n", n, it);
			#endif
			it++;
			nls.assemble();
			nls.solve(2, &Csln, &phisln);
			res_l2_norm = nls.get_residuum_l2_norm();

			#ifdef VERBAL
			info("Residuum L2 norm: %g\n", res_l2_norm);
			#endif

			Ci.copy(&Csln);
			phii.copy(&phisln);

		} while (res_l2_norm > NEWTON_TOL);
		#ifdef CONT_OUTPUT
		phiview.show(&phii);
		Cview.show(&Ci);
		#endif
		phip.copy(&phii);
		Cp.copy(&Ci);
	}
	verbose("\nTotal run time: %g sec", end_time());
	Cview.show(&Ci);
	//Cview.save_numbered_screenshot("screenshots/C%03d.bmp", REF_INIT, true);
	phiview.show(&phii);
	//phiview.save_numbered_screenshot("screenshots/phi%03d.bmp", REF_INIT, true);
	MeshView mview("small.mesh", 100, 30, 800, 800);
	mview.show(&mesh);
	//mview.save_numbered_screenshot("screenshots/mesh_refinelevel%03d.bmp", REF_INIT, true);
	View::wait();
}

void solveAdaptive(Mesh &Cmesh, Mesh &phimesh, Mesh &basemesh, NonlinSystem &nls, H1Space &C, H1Space &phi,
		Solution &Cp, Solution &Ci, Solution &phip, Solution &phii) {

	char title[100];
	//VectorView vview("electric field [V/m]", 0, 0, 600, 600);
	ScalarView Cview("Concentration [mol/m3]", 0, 0, 800, 800);
	ScalarView phiview("Voltage [V]", 650, 0, 600, 600);
	OrderView Cordview("C order", 0, 300, 600, 600);
	OrderView phiordview("Phi order", 600, 300, 600, 600);
	
	// Different Gnuplot graphs.

	// convergence graph wrt. the number of degrees of freedom
	GnuplotGraph graph_err;
	graph_err.set_captions("", "Timestep", "Error (Energy Norm)");
	graph_err.add_row("Reference solution", "k", "-", "O");

	// convergence graph wrt. CPU time
	GnuplotGraph graph_dof;
	graph_dof.set_captions("", "Timestep", "DOF");
	graph_dof.add_row(MULTIMESH ? "multi-mesh" : "single-mesh", "k", "-", "o");


	sprintf(title, "Initial iteration");
	phiview.set_title(title);
	Cview.set_title(title);
	phiview.show(&phii);
	Cview.show(&Ci);
	
	//Cview.save_numbered_screenshot("screenshots/C%03d.bmp", 0, true);
	//phiview.save_numbered_screenshot("screenshots/phi%03d.bmp", 0, true);
	
	Cordview.show(&C);
	phiordview.show(&phi);
	Solution Csln_coarse, phisln_coarse, Csln_fine, phisln_fine;
	int at_index = 1;	//for saving screenshot
	for (int n = 1; n <= NSTEP; n++) {
		if (n % UNREF_FREQ == 0) {
			Cmesh.copy(&basemesh);
			if (MULTIMESH) {
				phimesh.copy(&basemesh);
			}
			C.set_uniform_order(P_INIT);
			phi.set_uniform_order(P_INIT);
			int ndofs;
			ndofs = C.assign_dofs();
			phi.assign_dofs(ndofs);
		}	

		#ifdef VERBAL
		info("\n---- Time step %d ----", n);
		#endif

		int at = 0;
		bool done = false;
		double err;
		do {
			at++;
			at_index++;
			int it = 1;
			double res_l2_norm;
			if (n > 1 || at > 1) {
				nls.set_ic(&Csln_fine, &Ci, &phisln_fine, &phii);
			} else {
				/* No need to set anything, already set. */
				nls.set_ic(&Ci, &phii, &Ci, &phii);
			}
			//Loop for coarse mesh solution
			do {
				it++;

				#ifdef VERBAL
				info("\n -------- Time step %d, Newton iter %d --------\n", n, it);
				#endif

				nls.assemble();
				nls.solve(2, &Csln_coarse, &phisln_coarse);
				res_l2_norm = nls.get_residuum_l2_norm();

				#ifdef VERBAL
				info("Residuum L2 norm: %g\n", res_l2_norm);
				#endif

				Ci.copy(&Csln_coarse);
				phii.copy(&phisln_coarse);
			} while (res_l2_norm > NEWTON_TOL_COARSE);

			int ndf = C.get_num_dofs() + phi.get_num_dofs();
		    sprintf(title, "phi after COARSE solution, at=%d and n=%d, ndofs=%d", at, n, ndf);
		    phiview.set_title(title);
			phiview.show(&phii);
			sprintf(title, "C after COARSE solution, at=%d and n=%d, ndofs=%d. PRESS KEY TO CONTINUE", at, n, ndf);
			Cview.set_title(title);
			Cview.show(&Ci);
			Cview.wait_for_keypress();

			it = 1;
			// Loop for fine mesh solution
			RefNonlinSystem rs(&nls);
			rs.prepare();
			if (n > 1 || at > 1) {
				rs.set_ic(&Csln_fine, &phisln_fine, &Ci, &phii);
			} else {
				rs.set_ic(&Ci, &phii, &Ci, &phii);
			}

			do {
				it++;

				#ifdef VERBAL
				info("\n -------- Time step %d, Adaptivity step %d, Newton iter %d (Fine mesh) --------\n", n, at, it);
				#endif

				rs.assemble();
				rs.solve(2, &Csln_fine, &phisln_fine);
				res_l2_norm = rs.get_residuum_l2_norm();

				#ifdef VERBAL
				info("Residuum L2 norm: %g\n", res_l2_norm);
				#endif

				Ci.copy(&Csln_fine);
				phii.copy(&phisln_fine);
			} while (res_l2_norm > NEWTON_TOL_REF);

			ndf = C.get_num_dofs() + phi.get_num_dofs();
		    sprintf(title, "phi after FINE solution, at=%d and n=%d, ndofs=%d", at, n, ndf);
		    phiview.set_title(title);
			phiview.show(&phii);
			sprintf(title, "C after FINE solution, at=%d and n=%d, ndofs=%d. PRESS KEY TO CONTINUE", at, n, ndf);
			Cview.set_title(title);
			Cview.show(&Ci);
			Cview.wait_for_keypress();

			// Calculate element errors and total estimate
			H1OrthoHP hp(2, &C, &phi);
			info("\n Calculating element errors\n");
			err = hp.calc_error_2(&Csln_coarse, &phisln_coarse, &Csln_fine, &phisln_fine) * 100;
			info("Error: %g", err);

			if (err < ERR_STOP) {
				done = true;
			} else {
				hp.adapt(THRESHOLD, STRATEGY, H_ONLY);
			}

			int ndofs;
			ndofs = C.assign_dofs();
			ndofs += phi.assign_dofs(ndofs);
			info("NDOFS after adapting: %d", ndofs);
			if (ndofs >= MAX_NDOFS) {
				done = true;
			}

		    sprintf(title, "hp-mesh (C), time level %d, adaptivity %d", n, at);
		    Cordview.set_title(title);
		    sprintf(title, "hp-mesh (phi), time level %d, adaptivity %d", n, at);
		    phiordview.set_title(title);
		    Cordview.show(&C);
		    //Cordview.save_numbered_screenshot("screenshots/Cord%03d.bmp", at_index, true);
		    phiordview.show(&phi);
		    //phiordview.save_numbered_screenshot("screenshots/phiord%03d.bmp", at_index, true);


		} while (!done);
		phip.copy(&phii);
		Cp.copy(&Ci);

	   	graph_err.add_values(0, n, err);
		graph_err.save("error.gp");
		graph_dof.add_values(0, n, C.get_num_dofs() + phi.get_num_dofs());
		graph_dof.save("dofs.gp");
	    sprintf(title, "phi after time step %d", n);
	    phiview.set_title(title);
		phiview.show(&phii);
		sprintf(title, "C after time step %d", n);
		Cview.set_title(title);
		Cview.show(&Ci);
	}
	View::wait();
}

int main (int argc, char* argv[]) {

	bool adaptive = false;
	if (argc > 1) {
		std::string arg0(argv[1]);
		if (USE_ADAPTIVE.compare(arg0) == 0) {
			adaptive = true;
			info("Using adaptive solution");
		} else {
			error("Illegal argument %s", argv[1]);
			return -1;
		}
	} else {
		info("Using NON-adaptive solution");
	}

	// load the mesh file
	Mesh Cmesh, phimesh, basemesh;

	basemesh.load("small.mesh");
	basemesh.refine_towards_boundary(TOP_MARKER, REF_INIT);
	basemesh.refine_towards_boundary(BOT_MARKER, REF_INIT - 1);
	Cmesh.copy(&basemesh);
	phimesh.copy(&basemesh);

	// create the shapeset
	H1Shapeset shapeset;
	PrecalcShapeset Cpss(&shapeset);
	PrecalcShapeset phipss(&shapeset);

	// Spaces for concentration and the voltage
	H1Space C(&Cmesh, &shapeset);
	H1Space phi(MULTIMESH ? &phimesh : &Cmesh, &shapeset);

	// Initialize boundary conditions
	C.set_bc_types(C_bc_types);
	phi.set_bc_types(phi_bc_types);
	phi.set_bc_values(phi_bc_values);
	//C.set_bc_values(C_bc_values);
	
	// set polynomial degrees
	C.set_uniform_order(P_INIT);
	phi.set_uniform_order(P_INIT);
	
	
	// assign degrees of freedom
	int ndofs = 0;
	ndofs += C.assign_dofs(ndofs);
	ndofs += phi.assign_dofs(ndofs);
	info("ndofs: %d", ndofs);

	// The weak form for 2 equations
	WeakForm wf(2);

	Solution Cp,		// prveious time step solution, for the time integration
		Ci,		// solution convergin during the Newton's iteration
		phip,
		phii;

	
	// Add the bilinear and linear forms
	// generally, the equation system is described:
	// a11(u1, v1) + a12(u2, v1) + a1n(un, v1) = l1(v1)
	// a21(u1, v2) + a22(u2, v2) + a2n(un, v2) = l2(v2)
	// an1(u1, vn) + an2(u2, vn) + ann(un, vn) = ln(vn)
	wf.add_biform(0, 0, callback(J_euler_DFcDYc), UNSYM, ANY, 1, &phii);
	wf.add_biform(1, 1, callback(J_euler_DFphiDYphi), UNSYM);
	wf.add_biform(0, 1, callback(J_euler_DFcDYphi), UNSYM, ANY, 1, &Ci);
	wf.add_biform(1, 0, callback(J_euler_DFphiDYc), UNSYM);

	wf.add_liform(0, callback(Fc_euler), ANY, 3, &Cp, &Ci, &phii);
	wf.add_liform(1, callback(Fphi_euler), ANY, 2, &Ci, &phii);

	wf.add_liform_surf(1, callback(linear_form_surf_top), TOP_MARKER);
	
	// Noninear solver
	UmfpackSolver umfpack;
	NonlinSystem nls(&wf, &umfpack);
	nls.set_spaces(2, &C, &phi);
	if (MULTIMESH) {
		nls.set_pss(2, &Cpss, &phipss);
	} else {
		nls.set_pss(1, &Cpss);
	}
	
	info("UmfpackSolver initialized");

	// View initial guess for Newton's method
	// initial BC

	//Cp.set_dirichlet_lift(&C, &Cpss);
	//phip.set_dirichlet_lift(&phi, MULTIMESH ? &phipss : &Cpss);
	Cp.set_const(&Cmesh, C_CONC);
	phip.set_const(MULTIMESH ? &phimesh : &Cmesh, 0);

	Ci.copy(&Cp);
	phii.copy(&phip);

	nls.set_ic(&Ci, &phii, &Ci, &phii);

	if (adaptive) {
		solveAdaptive(Cmesh, phimesh, basemesh, nls, C, phi, Cp, Ci, phip, phii);
	} else {
		solveNonadaptive(Cmesh, nls, Cp, Ci, phip, phii);
	}
#define ERROR_SUCCESS 0
#define ERROR_FAILURE -1

	return ERROR_SUCCESS;
}
