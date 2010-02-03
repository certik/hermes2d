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

bool solveNonadaptive(Mesh &mesh, NonlinSystem &nls,
		Solution &Cp, Solution &Ci, Solution &phip, Solution &phii) {
	Solution Csln, phisln;
	for (int n = 1; n <= NSTEP; n++) {

		int it = 1;
		double res_l2_norm;

		do {

			it++;
			nls.assemble();
			nls.solve(2, &Csln, &phisln);
			res_l2_norm = nls.get_residuum_l2_norm();

			Ci.copy(&Csln);
			phii.copy(&phisln);

		} while (res_l2_norm > NEWTON_TOL);
		phip.copy(&phii);
		Cp.copy(&Ci);
	}
	scalar *sol_vector;
	int n_dof;
	nls.get_solution_vector(sol_vector, n_dof);
	printf("n_dof: %d\n", n_dof);
	double sum = 0;
	for (int i = 0; i < n_dof; i++) {
		 sum += sol_vector[i];
	}
	printf("coefficient sum = %g\n", sum);

	// Actual test. The value of 'sum' depend on the 
	// current shapeset. If you change the shapeset, 
	// you need to correct this number. 
	printf("ret: %g\n", fabs(sum - 3.90395e6));
	return !(fabs(sum - 3.90395e6) > 1);

}

void solveAdaptive(Mesh &Cmesh, Mesh &phimesh, Mesh &basemesh, NonlinSystem &nls, H1Space &C, H1Space &phi,
		Solution &Cp, Solution &Ci, Solution &phip, Solution &phii) {


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

		} while (!done);
		phip.copy(&phii);
		Cp.copy(&Ci);
	}
}

int main (int argc, char* argv[]) {
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

	//solveAdaptive(Cmesh, phimesh, basemesh, nls, C, phi, Cp, Ci, phip, phii);
	bool success = solveNonadaptive(Cmesh, nls, Cp, Ci, phip, phii);
	
	if (success) {
		printf("SUCCESSFUL\n");
	} else {
		printf("FAIL\n");
	}
	#define ERROR_SUCCESS 0
	#define ERROR_FAILURE -1

	return success ? ERROR_SUCCESS : ERROR_FAILURE;
}
