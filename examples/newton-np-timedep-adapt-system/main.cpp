#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"
#include <string>

// This example shows how to combine the automatic adaptivity with the
// Newton's method for a nonlinear time-dependent PDE system.
// The time discretization is done using implicit Euler or 
// Crank Nicholson method (see parameter TIME_DISCR).
// The following PDE's are solved:
// Nernst-Planck (describes the diffusion and migration of charged particles):
// dC/dt - D*div[grad(C)] - K*C*div[grad(phi)]=0
// where D and K are constants and C is the cation concentration variable,
// phi is the voltage variable in the Poisson equation:
// - div[grad(phi)] = L*(C - C0),
// where C0, and L are constant (anion concentration). C0 is constant
// anion concentration in the domain and L is material parameter.
// So, the equation variables are phi and C and the system describes the
// migration/diffusion of charged particles due to applied voltage.
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
// For Poisson equation, boundary 1 has a natural boundary condition 
// (electric field derivative is 0). 
// The voltage is applied to the boundaries 2 and 3 (Dirichlet boundaries)
// It is possible to adjust system paramter VOLT_BOUNDARY to apply
// Neumann boundary condition to 2 (instead of Dirichlet). But by default:
// BC 2: phi = VOLTAGE
// BC 3: phi = 0
// BC 1: dphi/dn = 0


#define SIDE_MARKER 1
#define TOP_MARKER 2
#define BOT_MARKER 3

// Parameters to tweak the amount of output to the console
#define NOSCREENSHOT

/*** Fundamental coefficients ***/
const double D = 10e-11; 	            // [m^2/s] Diffusion coefficient
const double R = 8.31; 		            // [J/mol*K] Gas constant
const double T = 293; 		            // [K] Aboslute temperature
const double F = 96485.3415;	        // [s * A / mol] Faraday constant
const double eps = 2.5e-2; 	          // [F/m] Electric permeability
const double mu = D / (R * T);        // Mobility of ions
const double z = 1;		                // Charge number
const double K = z * mu * F;          // Constant for equation
const double L =  F / eps;	          // Constant for equation
const double VOLTAGE = 1;	            // [V] Applied voltage
const scalar C0 = 1200;	              // [mol/m^3] Anion and counterion concentration

/* For Neumann boundary */
const double height = 180e-6;	              // [m] thickness of the domain
const double E_FIELD = VOLTAGE / height;    // Boundary condtion for positive voltage electrode


/* Simulation parameters */
const int PROJ_TYPE = 1;              // For the projection of the initial condition 
                                      // on the initial mesh: 1 = H1 projection, 0 = L2 projection
const int NSTEP = 50;                 // Number of time steps
const double TAU = 0.1;               // Size of the time step
const int P_INIT = 3;       	        // Initial polynomial degree of all mesh elements.
const int REF_INIT = 1;     	        // Number of initial refinements
const bool MULTIMESH = false;	        // Multimesh?
const int TIME_DISCR = 2;             // 1 for implicit Euler, 2 for Crank-Nicolson
const int VOLT_BOUNDARY = 1;          // 1 for Dirichlet, 2 for Neumann

/* Nonadaptive solution parameters */
const double NEWTON_TOL = 1e-6;       // Stopping criterion for nonadaptive solution

/* Adaptive solution parameters */
const double NEWTON_TOL_COARSE = 0.01;// Stopping criterion for Newton on coarse mesh
const double NEWTON_TOL_FINE = 0.05;  // Stopping criterion for Newton on fine mesh
const int NEWTON_MAX_ITER = 20;       // Maximum allowed number of Newton iterations
const bool NEWTON_ON_COARSE_MESH = false;  // true... Newton is done on coarse mesh in every adaptivity step
                                      // false...Newton is done on coarse mesh only once, then projection 
                                      // of the fine mesh solution to coarse mesh is used

const int UNREF_FREQ = 5;             // every UNREF_FREQth time step the mesh
                                      // is unrefined
const double THRESHOLD = 0.3;         // This is a quantitative parameter of the adapt(...) function and
                                      // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;               // Adaptive strategy:
                                      // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                      //   error is processed. If more elements have similar errors, refine
                                      //   all to keep the mesh symmetric.
                                      // STRATEGY = 1 ... refine all elements whose error is larger
                                      //   than THRESHOLD times maximum element error.
                                      // STRATEGY = 2 ... refine all elements whose error is larger
                                      //   than THRESHOLD.
                                      // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int ADAPT_TYPE = 0;             // Type of automatic adaptivity:
                                      // ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                      // ADAPT_TYPE = 1 ... adaptive h-FEM,
                                      // ADAPT_TYPE = 2 ... adaptive p-FEM.

const int NDOF_STOP = 5000;		        // To prevent adaptivity from going on forever.
const double ERR_STOP = 0.1;          // Stopping criterion for adaptivity (rel. error tolerance between the
                                      // fine mesh and coarse mesh solution in percent).

// Program parameters
const std::string USE_ADAPTIVE("adapt");

// Weak forms
#include "forms.cpp"


/*** Boundary types and conditions ***/

// Poisson takes Dirichlet and Neumann boundaries
int phi_bc_types(int marker) {
  return (marker == SIDE_MARKER || (marker == TOP_MARKER && VOLT_BOUNDARY == 2)) 
    ? BC_NATURAL : BC_ESSENTIAL;
}

// Nernst-Planck takes Neumann boundaries
int C_bc_types(int marker) {
  return BC_NATURAL;
}

// Diricleht Boundary conditions for Poisson equation.
scalar phi_bc_values(int marker, double x, double y) {
  return marker == TOP_MARKER ? VOLTAGE : 0.0;
}

template<class Real, class Scalar>
Scalar linear_form_surf_top(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
  return -E_FIELD * int_v<Real, Scalar>(n, wt, v);
}

scalar voltage_ic(double x, double y, double &dx, double &dy) {
  // y^2 function for the domain.
  //return (y+100e-6) * (y+100e-6) / (40000e-12);
  return 0.0;
}

scalar concentration_ic(double x, double y, double &dx, double &dy) {
  return C0;
}


/** Nonadaptive solver.*/
void solveNonadaptive(Mesh &mesh, NonlinSystem &nls,
    Solution &C_prev_time, Solution &C_prev_newton,
    Solution &phi_prev_time, Solution &phi_prev_newton) {

  //VectorView vview("electric field [V/m]", 0, 0, 600, 600);
  ScalarView Cview("Concentration [mol/m3]", 0, 0, 800, 800);
  ScalarView phiview("Voltage [V]", 650, 0, 600, 600);
  phiview.show(&phi_prev_time);
  Cview.show(&C_prev_time);
  char title[100];

  for (int n = 1; n <= NSTEP; n++) {
    verbose("\n---- Time step %d ----", n);
    if (!nls.solve_newton_2(&C_prev_newton, &phi_prev_newton,
         NEWTON_TOL, NEWTON_MAX_ITER)) error("Newton's method did not converge.");

    sprintf(title, "Solution, timestep = %i", n);
    phiview.set_title(title);
    phiview.show(&phi_prev_newton);
    Cview.set_title(title);
    Cview.show(&C_prev_newton);
    phi_prev_time.copy(&phi_prev_newton);
    C_prev_time.copy(&C_prev_newton);
    info("Just for debugging");
  }
  View::wait();
}

/** Adaptive solver.*/
void solveAdaptive(Mesh &Cmesh, Mesh &phimesh, Mesh &basemesh, NonlinSystem &nls,
     H1Space &Cspace, H1Space &phispace, Solution &C_prev_time, Solution &C_prev_newton,
     Solution &phi_prev_time, Solution &phi_prev_newton) {

  char title[100];
  //VectorView vview("electric field [V/m]", 0, 0, 600, 600);
  ScalarView Cview("Concentration [mol/m3]", 0, 0, 800, 800);
  ScalarView phiview("Voltage [V]", 650, 0, 600, 600);
  OrderView Cordview("C order", 0, 300, 600, 600);
  OrderView phiordview("Phi order", 600, 300, 600, 600);

  // convergence graph wrt. the number of degrees of freedom
  GnuplotGraph graph_err;
  graph_err.set_captions("", "Timestep", "Error (Energy Norm)");
  graph_err.add_row("Reference solution", "k", "-", "O");

  // convergence graph wrt. CPU time
  GnuplotGraph graph_dof;
  graph_dof.set_captions("", "Timestep", "DOF");
  graph_dof.add_row(MULTIMESH ? "multi-mesh" : "single-mesh", "k", "-", "o");


  phiview.set_title(title);
  Cview.set_title(title);
  phiview.show(&phi_prev_newton);
  Cview.show(&C_prev_newton);
  
  Cordview.show(&Cspace);
  phiordview.show(&phispace);
  
  //Newton's loop on a coarse mesh
  info("---- Time step 1, Newton solve on the coarse mesh\n");
  if (!nls.solve_newton_2(&C_prev_newton, &phi_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER)) 
    error("Newton's method did not converge.");

  Solution Csln_coarse, phisln_coarse, Csln_fine, phisln_fine;
  
  Csln_coarse.copy(&C_prev_newton);
  phisln_coarse.copy(&phi_prev_newton);

  int at_index = 1; //for saving screenshot
  int ndof;
  for (int n = 1; n <= NSTEP; n++) {
    if (n > 1 && n % UNREF_FREQ == 0) {
      Cmesh.copy(&basemesh);
      if (MULTIMESH) {
        phimesh.copy(&basemesh);
      }
      Cspace.set_uniform_order(P_INIT);
      phispace.set_uniform_order(P_INIT);
      ndof = assign_dofs(2, &Cspace, &phispace);

      // project the fine mesh solution on the globally derefined mesh
      info("---- Time step %d, projecting fine mesh solution on globally derefined mesh:\n", n);
      nls.set_ic(&Csln_fine, &phisln_fine, &C_prev_newton, &phi_prev_newton, PROJ_TYPE);

      if (NEWTON_ON_COARSE_MESH) {
        // Newton's loop on the globally derefined mesh
        info("---- Time step %d, Newton solve on globally derefined mesh:\n", n);
        if (!nls.solve_newton_2(&C_prev_newton, &phi_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER))
          error("Newton's method did not converge.");
      }
      Csln_coarse.copy(&C_prev_newton);
      phisln_coarse.copy(&phi_prev_newton);
    } 

    int at = 1;
    bool done = false;
    double err;
    do {
      // Loop for fine mesh solution
      RefNonlinSystem rs(&nls);
      rs.prepare();

      if (at == 1) rs.set_ic(&Csln_coarse, &phisln_coarse, &C_prev_newton, &phi_prev_newton);
      else rs.set_ic(&Csln_fine, &phisln_fine, &C_prev_newton, &phi_prev_newton);
      
      rs.solve_newton_2(&C_prev_newton, &phi_prev_newton, NEWTON_TOL_FINE, NEWTON_MAX_ITER);
      Csln_fine.copy(&C_prev_newton);
      phisln_fine.copy(&phi_prev_newton);

      sprintf(title, "phi Solution, time level %d, adapt step %d", n, at);
      phiview.set_title(title);
      AbsFilter mag(&phi_prev_newton);
      phiview.show(&mag); 
      sprintf(title, "C Solution, time level %d, adapt step %d", n, at);
      Cview.set_title(title);
      AbsFilter mag2(&C_prev_newton);
      Cview.show(&mag2); 

      // Calculate element errors and total estimate
      H1OrthoHP hp(2, &Cspace, &phispace);
      info("\n Calculating element errors\n");
      err = hp.calc_error_2(&Csln_coarse, &phisln_coarse, &Csln_fine, &phisln_fine) * 100;
      info("Error: %g%%", err);

      if (err < ERR_STOP) {
        done = true;
      } else {
        hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE);
        // enumerate degrees of freedom
        ndof = assign_dofs(2, &Cspace, &phispace);

        info("NDOF after adapting: %d", ndof);
        if (ndof >= NDOF_STOP) {
          info("NDOF reached to the max %d", NDOF_STOP);
          done = true;
        }

        // project the fine mesh solution on the new coarse mesh
        info("---- Time step %d, adaptivity step %d, projecting fine mesh solution on new coarse mesh:\n",
            n, at);
        nls.set_ic(&Csln_fine, &phisln_fine, &C_prev_newton, &phi_prev_newton, PROJ_TYPE);
        at++;
        if (NEWTON_ON_COARSE_MESH) {
          // Newton's loop on the globally derefined mesh
          info("---- Time step %d, Newton solve on globally derefined mesh:\n", n);
          if (!nls.solve_newton_2(&C_prev_newton, &phi_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER))
            error("Newton's method did not converge.");
        }

        Csln_coarse.copy(&C_prev_newton);
        phisln_coarse.copy(&phi_prev_newton);

      }
      
      sprintf(title, "hp-mesh (C), time level %d, adaptivity %d", n, at);
      Cordview.set_title(title);
      sprintf(title, "hp-mesh (phi), time level %d, adaptivity %d", n, at);
      phiordview.set_title(title);
      Cordview.show(&Cspace);
      phiordview.show(&phispace);
      #ifdef SCREENSHOT
      Cordview.save_numbered_screenshot("screenshots/Cord%03d.bmp", at_index, true);
      phiordview.save_numbered_screenshot("screenshots/phiord%03d.bmp", at_index, true);
      #endif
      at_index++;
    } while (!done);
    graph_err.add_values(0, n, err);
    graph_err.save("error.gp");
    graph_dof.add_values(0, n, Cspace.get_num_dofs() + phispace.get_num_dofs());
    graph_dof.save("dofs.gp");
    if (n == 1) {
      sprintf(title, "phi after time step %d, adjust the graph and PRESS ANY KEY", n);
    } else {
      sprintf(title, "phi after time step %d", n);
    }
    phiview.set_title(title);
    phiview.show(&phi_prev_newton);
    if (n == 1) {
      sprintf(title, "C after time step %d, adjust the graph and PRESS ANY KEY", n);
    } else {
      sprintf(title, "C after time step %d", n);
    }
    Cview.set_title(title);
    Cview.show(&C_prev_newton);
    #ifdef SCREENSHOT
    Cview.save_numbered_screenshot("screenshots/C%03d.bmp", n, true);
    phiview.save_numbered_screenshot("screenshots/phi%03d.bmp", n, true);
    #endif
    if (n == 1) {
      // Wait for key press, so one can go to 3D mode
      // which is way more informative in case of Nernst Planck
      View::wait(H2DV_WAIT_KEYPRESS);
    }
    phi_prev_time.copy(&phisln_fine);
    C_prev_time.copy(&Csln_fine);
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
  
  // When nonadaptive solution, refine the mesh
  basemesh.refine_towards_boundary(TOP_MARKER,
      adaptive ? REF_INIT : REF_INIT * 25);
  basemesh.refine_towards_boundary(BOT_MARKER,
    adaptive ? REF_INIT - 1 : (REF_INIT * 25) - 1);
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

  // set polynomial degrees
  C.set_uniform_order(P_INIT);
  phi.set_uniform_order(P_INIT);


  // enumerate degrees of freedom
  int ndof = assign_dofs(2, &C, &phi);

  info("ndof: %d", ndof);

  // The weak form for 2 equations
  WeakForm wf(2);

  Solution C_prev_time,    // prveious time step solution, for the time integration
    C_prev_newton,   // solution convergin during the Newton's iteration
    phi_prev_time,
    phi_prev_newton;

  // Add the bilinear and linear forms
  // generally, the equation system is described:
  if (TIME_DISCR == 1) {  //implicit euler
    wf.add_liform(0, callback(Fc_euler), H2D_ANY, 3,
        &C_prev_time, &C_prev_newton, &phi_prev_newton);
    wf.add_liform(1, callback(Fphi_euler), H2D_ANY, 2, &C_prev_newton, &phi_prev_newton);
    wf.add_biform(0, 0, callback(J_euler_DFcDYc), H2D_UNSYM, H2D_ANY, 1, &phi_prev_newton);
    wf.add_biform(0, 1, callback(J_euler_DFcDYphi), H2D_UNSYM, H2D_ANY, 1, &C_prev_newton);
    wf.add_biform(1, 0, callback(J_euler_DFphiDYc), H2D_UNSYM);
    wf.add_biform(1, 1, callback(J_euler_DFphiDYphi), H2D_UNSYM);
  } else {
    wf.add_liform(0, callback(Fc_cranic), H2D_ANY, 4,
        &C_prev_time, &C_prev_newton, &phi_prev_newton, &phi_prev_time);
    wf.add_liform(1, callback(Fphi_cranic), H2D_ANY, 2, &C_prev_newton, &phi_prev_newton);
    wf.add_biform(0, 0, callback(J_cranic_DFcDYc), H2D_UNSYM, H2D_ANY, 2, &phi_prev_newton, &phi_prev_time);
    wf.add_biform(0, 1, callback(J_cranic_DFcDYphi), H2D_UNSYM, H2D_ANY, 2, &C_prev_newton, &C_prev_time);
    wf.add_biform(1, 0, callback(J_cranic_DFphiDYc), H2D_UNSYM);
    wf.add_biform(1, 1, callback(J_cranic_DFphiDYphi), H2D_UNSYM);
    
  }
  // Neumann voltage boundary
  if (VOLT_BOUNDARY == 2) {
    wf.add_liform_surf(1, callback(linear_form_surf_top), TOP_MARKER);
  }

  // Nonlinear solver
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(2, &C, &phi);

  if (MULTIMESH) {
    nls.set_pss(2, &Cpss, &phipss);
  } else {
    nls.set_pss(1, &Cpss);
  }

  phi_prev_time.set_exact(MULTIMESH ? &phimesh : &Cmesh, voltage_ic);
  C_prev_time.set_exact(&Cmesh, concentration_ic);

  C_prev_newton.copy(&C_prev_time);
  phi_prev_newton.copy(&phi_prev_time);

  nls.set_ic(&C_prev_newton, &phi_prev_newton, &C_prev_newton, &phi_prev_newton);

  if (adaptive) {
    solveAdaptive(Cmesh, phimesh, basemesh, nls, C, phi, C_prev_time,
        C_prev_newton, phi_prev_time, phi_prev_newton);
  } else {
    solveNonadaptive(Cmesh, nls, C_prev_time, C_prev_newton, phi_prev_time, phi_prev_newton);
  }

  return 0;
}

