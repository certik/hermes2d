#include "hermes2d.h"

using namespace RefinementSelectors;

// This is a test for Nernst-Planck example

#define SIDE_MARKER 1
#define TOP_MARKER 2
#define BOT_MARKER 3

/*** Fundamental coefficients ***/
const double D = 1e-12; 	            // [m^2/s] Diffusion coefficient
const double R = 8.31; 		            // [J/mol*K] Gas constant
const double T = 293; 		            // [K] Aboslute temperature
const double F = 96485.3415;	            // [s * A / mol] Faraday constant
const double eps = 2.5e-2; 	            // [F/m] Electric permeability
const double mu = D / (R * T);              // Mobility of ions
const double z = 1;		            // Charge number
const double K = z * mu * F;                // Constant for equation
const double L =  F / eps;	            // Constant for equation
const double VOLTAGE = 1;	            // [V] Applied voltage
const double C_CONC = 1200;	            // [mol/m^3] Anion and counterion concentration

/* For Neumann boundary */
const double height = 180e-6;	            // [m] thickness of the domain
const double E_FIELD = VOLTAGE / height;    // Boundary condtion for positive voltage electrode


/* Simulation parameters */
const int NSTEP = 10;                // Number of time steps
const double TAU = 0.05;              // Size of the time step
const int P_INIT = 2;       	        // Initial polynomial degree of all mesh elements.
const int REF_INIT = 10;     	        // Number of initial refinements
const bool MULTIMESH = true;	        // Multimesh?

/* Nonadaptive solution parameters */
const double NEWTON_TOL = 1e-2;         // Stopping criterion for nonadaptive solution

/* Adaptive solution parameters */
const double NEWTON_TOL_COARSE = 0.01;  // Stopping criterion for Newton on coarse mesh
const double NEWTON_TOL_REF = 0.1;	    // Stopping criterion for Newton on fine mesh

const int UNREF_FREQ = 1;             // every UNREF_FREQth time step the mesh
                                        // is unrefined
const double THRESHOLD = 0.3;           // This is a quantitative parameter of the adapt(...) function and
                                        // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                 // Adaptive strategy:
                                        // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                        //   error is processed. If more elements have similar errors, refine
                                        //   all to keep the mesh symmetric.
                                        // STRATEGY = 1 ... refine all elements whose error is larger
                                        //   than THRESHOLD times maximum element error.
                                        // STRATEGY = 2 ... refine all elements whose error is larger
                                        //   than THRESHOLD.
                                        // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See the Sphinx tutorial (http://hpfem.org/hermes2d/doc/src/tutorial-2.html#adaptive-h-fem-and-hp-fem) for details.

const int NDOF_STOP = 5000;		          // To prevent adaptivity from going on forever.
const double ERR_STOP = 0.5;            // Stopping criterion for adaptivity (rel. error tolerance between the
                                        // fine mesh and coarse mesh solution in percent).

// Weak forms
#include "forms.cpp"


/*** Boundary types and conditions ***/

// Poisson takes Dirichlet and Neumann boundaries
BCType phi_bc_types(int marker) {
  return (marker == SIDE_MARKER || marker == TOP_MARKER)
    ? BC_NATURAL : BC_ESSENTIAL;
}

// Nernst-Planck takes Neumann boundaries
BCType C_bc_types(int marker) {
  return BC_NATURAL;
}

// Diricleht Boundary conditions for Poisson equation.
scalar essential_bc_values(int ess_bdy_marker, double x, double y) {
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
			nls.solve(Tuple<Solution*>(&Csln, &phisln));
			res_l2_norm = nls.get_residual_l2_norm();

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
	printf("ret: %g\n", fabs(sum - 50505));
	return !(fabs(sum - 50505) > 1);

}

bool solveAdaptive(Mesh &Cmesh, Mesh &phimesh, Mesh &basemesh, NonlinSystem &nls, H1Space &C, H1Space &phi,
                   Solution &Cp, Solution &Ci, Solution &phip, Solution &phii) {

  // create a selector which will select optimal candidate
  H1ProjBasedSelector selector(CAND_LIST, 1.0, H2DRS_DEFAULT_ORDER);

  Solution Csln_coarse, phisln_coarse, Csln_fine, phisln_fine;
  int at_index = 1; //for saving screenshot
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

    int at = 0;
    bool done = false;
    double err;
    do {
      at++;
      at_index++;
      int it = 1;
      double res_l2_norm;
      if (n > 1 || at > 1) {
        nls.project_global(Tuple<MeshFunction*>(&Csln_fine, &phisln_fine), 
                           Tuple<Solution*>(&Ci, &phii));
      } else {
        /* No need to set anything, already set. */
        nls.project_global(Tuple<MeshFunction*>(&Ci, &phii), 
                           Tuple<Solution*>(&Ci, &phii));
      }
      //Loop for coarse mesh solution
      do {
        it++;

        nls.assemble();
        nls.solve(Tuple<Solution*>(&Csln_coarse, &phisln_coarse));
        res_l2_norm = nls.get_residual_l2_norm();

        Ci.copy(&Csln_coarse);
        phii.copy(&phisln_coarse);
      } while (res_l2_norm > NEWTON_TOL_COARSE);

      it = 1;
      // Loop for fine mesh solution
      RefSystem rs(&nls);
      if (n > 1 || at > 1) {
        rs.project_global(Tuple<MeshFunction*>(&Csln_fine, &phisln_fine), 
                          Tuple<Solution*>(&Ci, &phii));
      } else {
        rs.project_global(Tuple<MeshFunction*>(&Ci, &phii), 
                          Tuple<Solution*>(&Ci, &phii));
      }

      do {
        it++;

        rs.assemble();
        rs.solve(Tuple<Solution*>(&Csln_fine, &phisln_fine));
        res_l2_norm = rs.get_residual_l2_norm();

        Ci.copy(&Csln_fine);
        phii.copy(&phisln_fine);
      } while (res_l2_norm > NEWTON_TOL_REF);

      // Calculate element errors and total estimate
      H1Adapt hp(&nls);
      hp.set_solutions(Tuple<Solution*>(&Csln_coarse, &phisln_coarse), Tuple<Solution*>(&Csln_fine, &phisln_fine));
      info("Calculating element errors");
      err = hp.calc_error() * 100;
      info("Error: %g%%", err);

      if (err < ERR_STOP) {
        done = true;
      } else {
        hp.adapt(&selector, THRESHOLD, STRATEGY);
      }

      int ndofs;
      ndofs = C.assign_dofs();
      ndofs += phi.assign_dofs(ndofs);
      info("NDOFS after adapting: %d", ndofs);
      if (ndofs >= NDOF_STOP) {
        info("NDOFs reached to the max %d", NDOF_STOP);
        done = true;
      }
    } while (!done);
    phip.copy(&phii);
    Cp.copy(&Ci);
  }
  int nd = C.get_num_dofs() + phi.get_num_dofs();
  info("NDOFs at the end of timestep: %d", nd);
	bool success = false;
  if (nd < 283)
    success = true;

	scalar *sol_vector;
  int n_dof;
	nls.get_solution_vector(sol_vector, n_dof);
	printf("n_dof: %d\n", n_dof);
	double sum = 0;
	for (int i = 0; i < n_dof; i++) {
		 sum += sol_vector[i];
	}
	printf("coefficient sum = %g\n", sum);

  return success;
}

int main (int argc, char* argv[]) {
  
  // load the mesh file
  Mesh Cmesh, phimesh, basemesh;

  H2DReader mloader;
  mloader.load("small.mesh", &basemesh);
  
  basemesh.refine_towards_boundary(TOP_MARKER, REF_INIT);
  basemesh.refine_towards_boundary(BOT_MARKER, REF_INIT - 1);
  Cmesh.copy(&basemesh);
  phimesh.copy(&basemesh);

  // Spaces for concentration and the voltage.
  H1Space C(&Cmesh, C_bc_types, NULL, P_INIT);
  H1Space phi(MULTIMESH ? &phimesh : &Cmesh, phi_bc_types, essential_bc_values, P_INIT);

  // set polynomial degrees.
  C.set_uniform_order(P_INIT);
  phi.set_uniform_order(P_INIT);
  info("ndof: %d", C.get_num_dofs() + phi.get_num_dofs());

  // The weak form for 2 equations.
  WeakForm wf(2);

  Solution C_prev_time,    // prveious time step solution, for the time integration
           C_prev_newton,   // solution convergin during the Newton's iteration
           phi_prev_time,
           phi_prev_newton;

  // Add the bilinear and linear forms
  // generally, the equation system is described:
  // a11(u1, v1) + a12(u2, v1) + a1n(un, v1) = l1(v1)
  // a21(u1, v2) + a22(u2, v2) + a2n(un, v2) = l2(v2)
  // an1(u1, vn) + an2(u2, vn) + ann(un, vn) = ln(vn)
  wf.add_matrix_form(0, 0, callback(J_euler_DFcDYc), H2D_UNSYM, H2D_ANY, &phi_prev_newton);
  wf.add_matrix_form(0, 1, callback(J_euler_DFcDYphi), H2D_UNSYM, H2D_ANY, &C_prev_newton);
  wf.add_matrix_form(1, 0, callback(J_euler_DFphiDYc), H2D_UNSYM);
  wf.add_matrix_form(1, 1, callback(J_euler_DFphiDYphi), H2D_UNSYM);

  wf.add_vector_form(0, callback(Fc_euler), H2D_ANY, Tuple<MeshFunction*>(&C_prev_time, &C_prev_newton, &phi_prev_newton));
  wf.add_vector_form(1, callback(Fphi_euler), H2D_ANY, Tuple<MeshFunction*>(&C_prev_newton, &phi_prev_newton));

  // Neumann voltage boundary.
  wf.add_vector_form_surf(1, callback(linear_form_surf_top), TOP_MARKER);

  // Nonlinear solver.
  NonlinSystem nls(&wf, Tuple<Space*>(&C, &phi));

  info("UmfpackSolver initialized");

  // View initial guess for Newton's method
  // initial BC
  
  C_prev_time.set_const(&Cmesh, C_CONC);
  phi_prev_time.set_const(MULTIMESH ? &phimesh : &Cmesh, 0);

  //  phi_prev_time.set_exact(MULTIMESH ? &phimesh : &Cmesh, voltage_ic);
  //  C_prev_time.set_exact(&Cmesh, concentration_ic);
  C_prev_newton.copy(&C_prev_time);
  phi_prev_newton.copy(&phi_prev_time);

  nls.project_global(Tuple<MeshFunction*>(&C_prev_newton, &phi_prev_newton), 
                     Tuple<Solution*>(&C_prev_newton, &phi_prev_newton));

  bool success = solveAdaptive(Cmesh, phimesh, basemesh, nls, C, phi, C_prev_time,
        C_prev_newton, phi_prev_time, phi_prev_newton);
  // bool success = solveNonadaptive(Cmesh, nls, C_prev_time, C_prev_newton, 
  // phi_prev_time, phi_prev_newton);

  if (success) {
    printf("SUCCESSFUL\n");
    } else {
    printf("FAIL\n");
  }
  #define ERROR_SUCCESS 0
  #define ERROR_FAILURE -1

  return success ? ERROR_SUCCESS : ERROR_FAILURE;
}
