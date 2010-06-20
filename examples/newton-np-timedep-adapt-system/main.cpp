#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"

#include "hermes2d.h"
#include <string>

using namespace RefinementSelectors;

/** \addtogroup e_newton_np_timedep_adapt_system Newton Time-dependant System with Adaptivity
 \{
 \brief This example shows how to combine the automatic adaptivity with the Newton's method for a nonlinear time-dependent PDE system.

 This example shows how to combine the automatic adaptivity with the
 Newton's method for a nonlinear time-dependent PDE system.
 The time discretization is done using implicit Euler or
 Crank Nicholson method (see parameter TIME_DISCR).
 The following PDE's are solved:
 Nernst-Planck (describes the diffusion and migration of charged particles):
 \f[dC/dt - D*div[grad(C)] - K*C*div[grad(\phi)]=0,\f]
 where D and K are constants and C is the cation concentration variable,
 phi is the voltage variable in the Poisson equation:
 \f[ - div[grad(\phi)] = L*(C - C_0),\f]
 where \f$C_0\f$, and L are constant (anion concentration). \f$C_0\f$ is constant
 anion concentration in the domain and L is material parameter.
 So, the equation variables are phi and C and the system describes the
 migration/diffusion of charged particles due to applied voltage.
 The simulation domain looks as follows:
 \verbatim
      2
  +----------+
  |          |
 1|          |1
  |          |
  +----------+
      3
 \endverbatim
 For the Nernst-Planck equation, all the boundaries are natural i.e. Neumann.
 Which basically means that the normal derivative is 0:
 \f[ BC: -D*dC/dn - K*C*d\phi/dn = 0 \f]
 For Poisson equation, boundary 1 has a natural boundary condition
 (electric field derivative is 0).
 The voltage is applied to the boundaries 2 and 3 (Dirichlet boundaries)
 It is possible to adjust system paramter VOLT_BOUNDARY to apply
 Neumann boundary condition to 2 (instead of Dirichlet). But by default:
  - BC 2: \f$\phi = VOLTAGE\f$
  - BC 3: \f$\phi = 0\f$
  - BC 1: \f$\frac{d\phi}{dn} = 0\f$
 */

#define SIDE_MARKER 1
#define TOP_MARKER 2
#define BOT_MARKER 3

// Parameters to tweak the amount of output to the console.
#define NOSCREENSHOT

/*** Fundamental coefficients ***/
const double D = 10e-11; 	            // [m^2/s] Diffusion coefficient.
const double R = 8.31; 		            // [J/mol*K] Gas constant.
const double T = 293; 		            // [K] Aboslute temperature.
const double F = 96485.3415;	            // [s * A / mol] Faraday constant.
const double eps = 2.5e-2; 	            // [F/m] Electric permeability.
const double mu = D / (R * T);              // Mobility of ions.
const double z = 1;		            // Charge number.
const double K = z * mu * F;                // Constant for equation.
const double L =  F / eps;	            // Constant for equation.
const double VOLTAGE = 1;	            // [V] Applied voltage.
const scalar C0 = 1200;	                    // [mol/m^3] Anion and counterion concentration.

/* For Neumann boundary */
const double height = 180e-6;	            // [m] thickness of the domain.
const double E_FIELD = VOLTAGE / height;    // Boundary condtion for positive voltage electrode


/* Simulation parameters */
const int NSTEP = 50;                 // Number of time steps.
const double TAU = 0.1;               // Size of the time step.
const int P_INIT = 3;       	      // Initial polynomial degree of all mesh elements.
const int REF_INIT = 1;     	      // Number of initial refinements.
const bool MULTIMESH = false;	      // Multimesh?
const int TIME_DISCR = 2;             // 1 for implicit Euler, 2 for Crank-Nicolson.
const int VOLT_BOUNDARY = 1;          // 1 for Dirichlet, 2 for Neumann.

/* Nonadaptive solution parameters */
const double NEWTON_TOL = 1e-6;       // Stopping criterion for nonadaptive solution.

/* Adaptive solution parameters */
const bool SOLVE_ON_COARSE_MESH = false;  // true... Newton is done on coarse mesh in every adaptivity step.
                                      // false...Newton is done on coarse mesh only once, then projection 
                                      // of the fine mesh solution to coarse mesh is used.
const double NEWTON_TOL_COARSE = 0.01;// Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.05;  // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 100;      // Maximum allowed number of Newton iterations.

const int UNREF_FREQ = 5;             // every UNREF_FREQth time step the mesh is unrefined.
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
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;  // Maximum allowed level of hanging nodes:
                                 // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                 // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                 // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                 // Note that regular meshes are not supported, this is due to
                                 // their notoriously bad performance.
const double CONV_EXP = 1.0;     // Default value is 1.0. This parameter influences the selection of
                                 // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int NDOF_STOP = 5000;		        // To prevent adaptivity from going on forever.
const double ERR_STOP = 0.1;          // Stopping criterion for adaptivity (rel. error tolerance between the
                                      // fine mesh and coarse mesh solution in percent).

// Program parameters
const std::string USE_ADAPTIVE("adapt");

// Weak forms
#include "forms.cpp"


/*** Boundary types and conditions ***/

// Poisson takes Dirichlet and Neumann boundaries
BCType phi_bc_types(int marker) {
  return (marker == SIDE_MARKER || (marker == TOP_MARKER && VOLT_BOUNDARY == 2)) 
    ? BC_NATURAL : BC_ESSENTIAL;
}

// Nernst-Planck takes Neumann boundaries
BCType C_bc_types(int marker) {
  return BC_NATURAL;
}

// Diricleht Boundary conditions for Poisson equation.
scalar C_essential_bc_values(int ess_bdy_marker, double x, double y) {
  return 0;
}

// Diricleht Boundary conditions for Poisson equation.
scalar phi_essential_bc_values(int ess_bdy_marker, double x, double y) {
  return ess_bdy_marker == TOP_MARKER ? VOLTAGE : 0.0;
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
    verbose("---- Time step %d:", n);
    bool verbose = true; // Default is false.
    if (!nls.solve_newton(Tuple<Solution*>(&C_prev_newton, &phi_prev_newton),
			  NEWTON_TOL, NEWTON_MAX_ITER, verbose)) 
      error("Newton's method did not converge.");

    sprintf(title, "Solution, timestep = %i", n);
    phiview.set_title(title);
    phiview.show(&phi_prev_newton);
    Cview.set_title(title);
    Cview.show(&C_prev_newton);
    phi_prev_time.copy(&phi_prev_newton);
    C_prev_time.copy(&C_prev_newton);
    //info("Just for debugging");
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

  // Convergence graph wrt. the number of degrees of freedom.
  SimpleGraph graph_err, graph_dof;

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Visualization.
  phiview.set_title(title);
  Cview.set_title(title);
  phiview.show(&phi_prev_newton);
  Cview.show(&C_prev_newton);  
  Cordview.show(&Cspace);
  phiordview.show(&phispace);
  
  //Newton's loop on a coarse mesh
  info("---- Time step 1, Newton solve on the coarse mesh\n");
  bool verbose = true; // Default is false.
  if (!nls.solve_newton(Tuple<Solution*>(&C_prev_newton, &phi_prev_newton), 
                        NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) 
    error("Newton's method did not converge.");

  Solution Csln_coarse, phisln_coarse, Csln_fine, phisln_fine;
  
  Csln_coarse.copy(&C_prev_newton);
  phisln_coarse.copy(&phi_prev_newton);

  int at_index = 1; // For saving screenshots.
  int ndof;
  for (int n = 1; n <= NSTEP; n++) {
    if (n > 1 && n % UNREF_FREQ == 0) {
      Cmesh.copy(&basemesh);
      if (MULTIMESH) {
        phimesh.copy(&basemesh);
      }
      Cspace.set_uniform_order(P_INIT);
      phispace.set_uniform_order(P_INIT);

      // Project the fine mesh solution on the globally derefined mesh.
      info("---- Time step %d, projecting fine mesh solution on globally derefined mesh:\n", n);
      nls.project_global(Tuple<MeshFunction*>(&Csln_fine, &phisln_fine), 
                         Tuple<Solution*>(&C_prev_newton, &phi_prev_newton));

      if (SOLVE_ON_COARSE_MESH) {
        // Newton's loop on the globally derefined mesh.
        info("---- Time step %d, Newton solve on globally derefined mesh:\n", n);
        if (!nls.solve_newton(Tuple<Solution*>(&C_prev_newton, &phi_prev_newton), 
                              NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
          error("Newton's method did not converge.");
      }
      Csln_coarse.copy(&C_prev_newton);
      phisln_coarse.copy(&phi_prev_newton);
    } 

    int as = 1;
    bool done = false;
    double err_est;
    do {
      // Loop for fine mesh solution.
      RefSystem rs(&nls);

      // Set initial condition for the Newton's method on the fine mesh.
      if (as == 1) {         
        info("Projecting coarse mesh solutions to obtain initial vector on new fine meshes.");
        rs.project_global(Tuple<MeshFunction*>(&Csln_coarse, &phisln_coarse), 
                          Tuple<Solution*>(&C_prev_newton, &phi_prev_newton));
      }
      else {
        info("Projecting previous fine mesh solutions to obtain initial vector on new fine meshes.");
        rs.project_global(Tuple<MeshFunction*>(&Csln_fine, &phisln_fine), 
                          Tuple<Solution*>(&C_prev_newton, &phi_prev_newton));
      }
      
      // Newton's method on fine mesh
      info("Solving on fine mesh.");
      rs.solve_newton(Tuple<Solution*>(&C_prev_newton, &phi_prev_newton), 
                      NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose);

      // Store the result in sln_fine.
      Csln_fine.copy(&C_prev_newton);
      phisln_fine.copy(&phi_prev_newton);

      // Visualization.
      sprintf(title, "phi Solution, time level %d, adapt step %d", n, as);
      phiview.set_title(title);
      AbsFilter mag(&phi_prev_newton);
      phiview.show(&mag); 
      sprintf(title, "C Solution, time level %d, adapt step %d", n, as);
      Cview.set_title(title);
      AbsFilter mag2(&C_prev_newton);
      Cview.show(&mag2); 

      // Calculate error estimate wrt. fine mesh solution.
      info("Calculating error.");
      H1Adapt hp(&nls);
      hp.set_solutions(Tuple<Solution*>(&Csln_coarse, &phisln_coarse), 
                       Tuple<Solution*>(&Csln_fine, &phisln_fine));
      err_est = hp.calc_error() * 100;
      info("Rel. error: %g%%", err_est);

      // If err_est too large, adapt the mesh.
      if (err_est < ERR_STOP) {
        done = true;
      } 
      else {
        hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
        info("NDOF after adapting: %d", nls.get_num_dofs());
        if (nls.get_num_dofs() >= NDOF_STOP) {
          info("NDOF reached to the max %d", NDOF_STOP);
          done = true;
        }

        // Project the fine mesh solution on the new coarse mesh.
        info("---- Time step %d, adaptivity step %d, projecting fine mesh solution on new coarse mesh:\n",
            n, as);

        // Project the fine mesh solution on the new coarse mesh.
        if (SOLVE_ON_COARSE_MESH) 
          info("Projecting fine mesh solutions to obtain initial vector on new coarse meshes.");
        else 
          info("Projecting fine mesh solutions on coarse meshes for error calculation.");
        nls.project_global(Tuple<MeshFunction*>(&Csln_fine, &phisln_fine), 
                           Tuple<Solution*>(&C_prev_newton, &phi_prev_newton));


        if (SOLVE_ON_COARSE_MESH) {
          // Newton's loop on the coarse mesh.
          info("---- Time step %d, solving on new coarse meshes:\n", n);
          if (!nls.solve_newton(Tuple<Solution*>(&C_prev_newton, &phi_prev_newton), 
                                NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
            error("Newton's method did not converge.");
        }

        // Store the result in sln_coarse.
        Csln_coarse.copy(&C_prev_newton);
        phisln_coarse.copy(&phi_prev_newton);

        as++;
      }
      
      // Visualize the solution and mesh.
      sprintf(title, "hp-mesh (C), time level %d, adaptivity %d", n, as);
      Cordview.set_title(title);
      sprintf(title, "hp-mesh (phi), time level %d, adaptivity %d", n, as);
      phiordview.set_title(title);
      Cordview.show(&Cspace);
      phiordview.show(&phispace);
      #ifdef SCREENSHOT
      Cordview.save_numbered_screenshot("screenshots/Cord%03d.bmp", at_index, true);
      phiordview.save_numbered_screenshot("screenshots/phiord%03d.bmp", at_index, true);
      #endif
      at_index++;
    } while (!done);

    graph_err.add_values(n, err_est);
    graph_err.save("error.gp");
    graph_dof.add_values(n, nls.get_num_dofs());
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
      // which is way more informative in case of Nernst Planck.
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

  // Load the mesh file.
  Mesh Cmesh, phimesh, basemesh;
  H2DReader mloader;
  mloader.load("small.mesh", &basemesh);
  
  // When nonadaptive solution, refine the mesh.
  basemesh.refine_towards_boundary(TOP_MARKER,
      adaptive ? REF_INIT : REF_INIT * 25);
  basemesh.refine_towards_boundary(BOT_MARKER,
    adaptive ? REF_INIT - 1 : (REF_INIT * 25) - 1);
  Cmesh.copy(&basemesh);
  phimesh.copy(&basemesh);

  // Spaces for concentration and the voltage.
  H1Space C(&Cmesh, C_bc_types, C_essential_bc_values, P_INIT);
  H1Space phi(MULTIMESH ? &phimesh : &Cmesh, phi_bc_types, phi_essential_bc_values, P_INIT);

  // The weak form for 2 equations.
  WeakForm wf(2);

  Solution C_prev_time,    // Previous time step solution, for the time integration.
    C_prev_newton,         // Solution converging during the Newton's iteration.
    phi_prev_time,
    phi_prev_newton;

  // Add the bilinear and linear forms.
  if (TIME_DISCR == 1) {  // Implicit Euler.
    wf.add_vector_form(0, callback(Fc_euler), H2D_ANY,
		  Tuple<MeshFunction*>(&C_prev_time, &C_prev_newton, &phi_prev_newton));
    wf.add_vector_form(1, callback(Fphi_euler), H2D_ANY, Tuple<MeshFunction*>(&C_prev_newton, &phi_prev_newton));
    wf.add_matrix_form(0, 0, callback(J_euler_DFcDYc), H2D_UNSYM, H2D_ANY, &phi_prev_newton);
    wf.add_matrix_form(0, 1, callback(J_euler_DFcDYphi), H2D_UNSYM, H2D_ANY, &C_prev_newton);
    wf.add_matrix_form(1, 0, callback(J_euler_DFphiDYc), H2D_UNSYM);
    wf.add_matrix_form(1, 1, callback(J_euler_DFphiDYphi), H2D_UNSYM);
  } else {
    wf.add_vector_form(0, callback(Fc_cranic), H2D_ANY, 
		  Tuple<MeshFunction*>(&C_prev_time, &C_prev_newton, &phi_prev_newton, &phi_prev_time));
    wf.add_vector_form(1, callback(Fphi_cranic), H2D_ANY, Tuple<MeshFunction*>(&C_prev_newton, &phi_prev_newton));
    wf.add_matrix_form(0, 0, callback(J_cranic_DFcDYc), H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&phi_prev_newton, &phi_prev_time));
    wf.add_matrix_form(0, 1, callback(J_cranic_DFcDYphi), H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&C_prev_newton, &C_prev_time));
    wf.add_matrix_form(1, 0, callback(J_cranic_DFphiDYc), H2D_UNSYM);
    wf.add_matrix_form(1, 1, callback(J_cranic_DFphiDYphi), H2D_UNSYM);
  }

  // Neumann voltage boundary.
  if (VOLT_BOUNDARY == 2) {
    wf.add_vector_form_surf(1, callback(linear_form_surf_top), TOP_MARKER);
  }

  NonlinSystem nls(&wf, Tuple<Space*>(&C, &phi));

  // Set initial conditions.
  phi_prev_time.set_exact(MULTIMESH ? &phimesh : &Cmesh, voltage_ic);
  C_prev_time.set_exact(&Cmesh, concentration_ic);

  // Project the function init_cond() on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  info("Projecting initial conditions to obtain initial vector for the Newton'w method.");
  nls.project_global(Tuple<MeshFunction*>(&C_prev_time, &phi_prev_time), 
                     Tuple<Solution*>(&C_prev_newton, &phi_prev_newton));

  if (adaptive) {
    solveAdaptive(Cmesh, phimesh, basemesh, nls, C, phi, C_prev_time,
        C_prev_newton, phi_prev_time, phi_prev_newton);
  } else {
    solveNonadaptive(Cmesh, nls, C_prev_time, C_prev_newton, phi_prev_time, phi_prev_newton);
  }

  return 0;
}
/// \}
