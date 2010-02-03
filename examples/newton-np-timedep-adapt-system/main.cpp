#include "hermes2d.h"
#include "solver_umfpack.h"
#include <string>

// This example shows how to combine the automatic adaptivity with the
// Newton's method for a nonlinear time-dependent PDE system.
// The time discretization is done using implicit Euler method.
// Some problem parameters can be changed below.
// The following PDE's are solved:
// Nernst-Planck (describes the diffusion and migration of charged particles):
// dC/dt - D*div[grad(C)] - K*C*div[grad(phi)]=0
// Poisson equation (describes the electric field):
// - div[grad(phi)] = L*(C - C0)
// The equation variables are phi and C and the system describes the
// migration/diffusion of charged particles due to electric field.
// The simulation domain looks as follows:
//      2
//  ____________
//  |          |
// 1|          |1
//  ____________
//      3
// For the Nernst-Planck equation, all the boundaries are natural i.e. Neumann.
// Which basically means that the normal derivative is 0:
// BC: -D*dC/dn - K*C*dphi/dn = 0
// For Poisson equation, 1 has natural boundary conditions (electric field
// derivative is 0). The voltage is applied to the boundaries 2 and 3.
// However, to make the equations stable, i.e. maintain the charge balance
// inside the domain, the positive voltages on the boundary 2 must be
// described by using Neumann boundary. Boundary 3 is an essential boundary, i.e:
// BC 2: dphi/dn = E_FIELD
// BC 3: phi = 0
// BC 1: dphi/dn = 0


#define SIDE_MARKER 1
#define TOP_MARKER 2
#define BOT_MARKER 3

// Parameters to tweak the amount of output to the console
#define VERBOSE
#define NONCONT_OUTPUT

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
const int NSTEP = 20;                   // Number of time steps
const double TAU = 0.05;                // Size of the time step
const int P_INIT = 2;       	        // Initial polynomial degree of all mesh elements.
const int REF_INIT = 12;     	        // Number of initial refinements
const bool MULTIMESH = true;	        // Multimesh?

/* Nonadaptive solution parameters */
const double NEWTON_TOL = 1e-2;         // Stopping criterion for nonadaptive solution

/* Adaptive solution parameters */
const double NEWTON_TOL_COARSE = 0.05;  // Stopping criterion for Newton on coarse mesh
const double NEWTON_TOL_REF = 0.5;	// Stopping criterion for Newton on fine mesh

const int UNREF_FREQ = 10;           	// every UNREF_FREQth time step the mesh
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
const int ADAPT_TYPE = 0;               // Type of automatic adaptivity:
                                        // ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                        // ADAPT_TYPE = 1 ... adaptive h-FEM,
                                        // ADAPT_TYPE = 2 ... adaptive p-FEM.

const int NDOF_STOP = 5000;		// To prevent adaptivity from going on forever.
const double ERR_STOP = 0.5;            // Stopping criterion for adaptivity (rel. error tolerance between the
                                        // fine mesh and coarse mesh solution in percent).

// Program parameters 
const std::string USE_ADAPTIVE("adapt");

// Weak forms
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


/** Nonadaptive solver.*/
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

    #ifdef VERBOSE
    info("\n---- Time step %d ----", n);
    #endif

    int it = 1;
    double res_l2_norm;

    do {

      #ifdef VERBOSE
      info("\n -------- Time step %d, Newton iter %d --------\n", n, it);
      #endif
      it++;
      nls.assemble();
      nls.solve(2, &Csln, &phisln);
      res_l2_norm = nls.get_residuum_l2_norm();

      #ifdef VERBOSE
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
  phiview.show(&phii);
  //MeshView mview("small.mesh", 100, 30, 800, 800);
  //mview.show(&mesh);
  View::wait();
}

/** Adaptive solver.*/
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
  
  Cordview.show(&C);
  phiordview.show(&phi);
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

    #ifdef VERBOSE
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

        #ifdef VERBOSE
        info("\n -------- Time step %d, Newton iter %d --------\n", n, it);
        #endif

        nls.assemble();
        nls.solve(2, &Csln_coarse, &phisln_coarse);
        res_l2_norm = nls.get_residuum_l2_norm();

        #ifdef VERBOSE
        info("Residuum L2 norm: %g\n", res_l2_norm);
        #endif

        Ci.copy(&Csln_coarse);
        phii.copy(&phisln_coarse);
      } while (res_l2_norm > NEWTON_TOL_COARSE);

      int ndf = C.get_num_dofs() + phi.get_num_dofs();
      sprintf(title, "phi after COARSE solution, at=%d and n=%d, ndofs=%d", at, n, ndf);
      phiview.set_title(title);
      phiview.show(&phii);
      sprintf(title, "C after COARSE solution, at=%d and n=%d,\
           ndofs=%d. PRESS KEY TO CONTINUE", at, n, ndf);
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

        #ifdef VERBOSE
        info("\n -------- Time step %d, Adaptivity step %d, Newton iter %d (Fine mesh) --------\n", n, at, it);
        #endif

        rs.assemble();
        rs.solve(2, &Csln_fine, &phisln_fine);
        res_l2_norm = rs.get_residuum_l2_norm();

        #ifdef VERBOSE
        info("Residuum L2 norm: %g\n", res_l2_norm);
        #endif

        Ci.copy(&Csln_fine);
        phii.copy(&phisln_fine);
      } while (res_l2_norm > NEWTON_TOL_REF);

      ndf = C.get_num_dofs() + phi.get_num_dofs();
      sprintf(title, "phi after FINE solution, at=%d and n=%d, ndofs=%d", at, n, ndf);
      phiview.set_title(title);
      phiview.show(&phii);
      sprintf(title, "C after FINE solution, at=%d and n=%d, ndofs=%d.\
           PRESS KEY TO CONTINUE", at, n, ndf);
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
        hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE);
      }

      int ndofs;
      ndofs = C.assign_dofs();
      ndofs += phi.assign_dofs(ndofs);
      info("NDOFS after adapting: %d", ndofs);
      if (ndofs >= NDOF_STOP) {
        done = true;
      }

      sprintf(title, "hp-mesh (C), time level %d, adaptivity %d", n, at);
      Cordview.set_title(title);
      sprintf(title, "hp-mesh (phi), time level %d, adaptivity %d", n, at);
      phiordview.set_title(title);
      Cordview.show(&C);
      phiordview.show(&phi);
      #ifdef SCREENSHOT
      Cordview.save_numbered_screenshot("screenshots/Cord%03d.bmp", at_index, true);
      phiordview.save_numbered_screenshot("screenshots/phiord%03d.bmp", at_index, true);
      #endif

      

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

  H2DReader mloader;
  mloader.load("small.mesh", &basemesh);
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

  Solution Cp,    // prveious time step solution, for the time integration
    Ci,   // solution convergin during the Newton's iteration
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
  
  // Nonlinear solver
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(2, &C, &phi);
  if (MULTIMESH) {
    nls.set_pss(2, &Cpss, &phipss);
  } else {
    nls.set_pss(1, &Cpss);
  }
  
  info("UmfpackSolver initialized");


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

  return 0;
}

