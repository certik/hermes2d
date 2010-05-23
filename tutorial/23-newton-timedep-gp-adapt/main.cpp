#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

using namespace RefinementSelectors;

//  This example shows how to combine automatic adaptivity with the Newton's
//  method for a nonlinear complex-valued time-dependent PDE (the Gross-Pitaevski
//  equation describing the behavior of Einstein-Bose quantum gases)
//  discretized implicitly in time (via implicit Euler or Crank-Nicolson).
//
//  PDE: non-stationary complex Gross-Pitaevski equation
//  describing resonances in Bose-Einstein condensates
//
//  ih \partial \psi/\partial t = -h^2/(2m) \Delta \psi +
//  g \psi |\psi|^2 + 1/2 m \omega^2 (x^2 + y^2) \psi
//
//  square (-1, 1)^2
//
//  BC:  homogeneous Dirichlet everywhere on the boundary
//
//  Time-stepping: either implicit Euler or Crank-Nicolson

const bool NEWTON_ON_COARSE_MESH = false;  // true... Newton is done on coarse mesh in every adaptivity step
                                           // false...Newton is done on coarse mesh only once, then projection
                                           // of the fine mesh solution to coarse mesh is used
const int P_INIT = 1;                      // Initial polynomial degree
const int PROJ_TYPE = 1;                   // For the projection of the initial condition
                                           // on the initial mesh: 1 = H1 projection, 0 = L2 projection
const int TIME_DISCR = 2;                  // 1 for implicit Euler, 2 for Crank-Nicolson
const double T_FINAL = 200.0;              // Time interval length
const double TAU = 0.01;                   // Time step
const int INIT_REF_NUM = 2;                // Number of initial uniform refinements

// Adaptivity
const int UNREF_FREQ = 1;                  // Every UNREF_FREQ time step the mesh is unrefined
const double THRESHOLD = 0.3;              // This is a quantitative parameter of the adapt(...) function and
                                           // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                    // Adaptive strategy:
                                           // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                           //   error is processed. If more elements have similar errors, refine
                                           //   all to keep the mesh symmetric.
                                           // STRATEGY = 1 ... refine all elements whose error is larger
                                           //   than THRESHOLD times maximum element error.
                                           // STRATEGY = 2 ... refine all elements whose error is larger
                                           //   than THRESHOLD.
                                           // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;   // Predefined list of element refinement candidates. Possible values are
                                           // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                           // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                           // See the User Documentation for details.
const int MESH_REGULARITY = -1;            // Maximum allowed level of hanging nodes:
                                           // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                           // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                           // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                           // Note that regular meshes are not supported, this is due to
                                           // their notoriously bad performance.
const double CONV_EXP = 1.0;               // Default value is 1.0. This parameter influences the selection of
                                           // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int MAX_ORDER = 5;                   // Maximum polynomial order allowed in hp-adaptivity
                                           // had to be limited due to complicated integrals
const double ERR_STOP = 2.0;               // Stopping criterion for hp-adaptivity
                                           // (relative error between reference and coarse solution in percent)
const int NDOF_STOP = 60000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.

// Newton's method
const double NEWTON_TOL_COARSE = 0.01;     // Stopping criterion for Newton on coarse mesh
const double NEWTON_TOL_FINE = 0.05;       // Stopping criterion for Newton on fine mesh
const int NEWTON_MAX_ITER = 20;            // Maximum allowed number of Newton iterations

// Problem constants
const double H = 1;                      // Planck constant 6.626068e-34;
const double M = 1;                      // mass of boson
const double G = 1;                      // coupling constant
const double OMEGA = 1;                  // frequency


// Initial conditions
scalar fn_init(double x, double y, scalar& dx, scalar& dy)
{
  scalar val = exp(-20*(x*x + y*y));
  dx = val * (-40.0*x);
  dy = val * (-40.0*y);
  return val;
}

// Boundary condition types
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Boundary condition values
scalar bc_values(int marker, double x, double y)
{
 return 0;
}

// Weak forms
# include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh file.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("square.mesh", &basemesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

  // Initialize the shapeset and the cache.
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // Create finite element space.
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // Enumerate degrees of freedom.
  int ndof = assign_dofs(&space);

  // Solutions for the Newton's iteration and adaptivity.
  Solution Psi_prev_time, Psi_prev_newton;

  // Initialize the weak formulation.
  WeakForm wf;
  if(TIME_DISCR == 1) {
    wf.add_biform(callback(jacobian_euler), H2D_UNSYM, H2D_ANY, 1, &Psi_prev_newton);
    wf.add_liform(callback(residual_euler), H2D_ANY, 2, &Psi_prev_newton, &Psi_prev_time);
  }
  else {
    wf.add_biform(callback(jacobian_cranic), H2D_UNSYM, H2D_ANY, 1, &Psi_prev_newton);
    wf.add_liform(callback(residual_cranic), H2D_ANY, 2, &Psi_prev_newton, &Psi_prev_time);
  }

  // Initialize the nonlinear system and solver.
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_space(&space);
  nls.set_pss(&pss);

  // Visualization.
  char title[100];
  ScalarView view("", 0, 0, 600, 500);
  view.fix_scale_width(80);
  ScalarView magview("", 0, 0, 600, 500);
  magview.fix_scale_width(80);
  OrderView ordview("", 610, 0, 600, 500);
  ordview.fix_scale_width(80);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_time_dof, graph_time_err;

  // Set initial condition at zero time level.
  Psi_prev_time.set_exact(&mesh, fn_init);

  // Project fn_init() on the FE space and use it as initial 
  // condition for the Newton's method.
  nls.project_global(&fn_init, &Psi_prev_newton, PROJ_TYPE);

  // Newton's loop on the coarse mesh.
  info("---- Time step 1, Newton solve on coarse mesh:");
  if (!nls.solve_newton(&Psi_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER))
    error("Newton's method did not converge.");

  // Store the result in sln_coarse.
  Solution sln_coarse, sln_fine;
  sln_coarse.copy(&Psi_prev_newton);

  // Create selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, MAX_ORDER, &shapeset);

  // Time stepping loop.
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int n = 1; n <= nstep; n++)
  {
    // Periodic global derefinements.
    if (n > 1 && n % UNREF_FREQ == 0) {
      info("---- Time step %d, global derefinement.", n);
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);
      ndof = assign_dofs(&space);

      // Project the fine mesh solution on the globally derefined mesh.
      info("---- Time step %d, projecting fine mesh solution on globally derefined mesh:", n);
      nls.project_global(&sln_fine, &Psi_prev_newton, PROJ_TYPE);

      if (NEWTON_ON_COARSE_MESH) {
        // Newton's loop on the globally derefined mesh.
        info("---- Time step %d, Newton solve on globally derefined mesh:", n);
        if (!nls.solve_newton(&Psi_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER))
          error("Newton's method did not converge.");
      }

      // Store the result in sln_coarse.
      sln_coarse.copy(&Psi_prev_newton);
    }

    // Adaptivity loop.
    bool done = false;
    double err_est;
    int as = 1;
    do
    {
      info("---- Time step %d, adaptivity step %d, Newton solve on fine mesh:", n, as);

      // Reference nonlinear system.
      RefNonlinSystem rnls(&nls);
      rnls.prepare();

      // Set initial condition for the Newton's method on the fine mesh.
      if (as == 1) rnls.project_global(&sln_coarse, &Psi_prev_newton);
      else rnls.project_global(&sln_fine, &Psi_prev_newton);

      // Newton's method on fine mesh.
      if (!rnls.solve_newton(&Psi_prev_newton, NEWTON_TOL_FINE, NEWTON_MAX_ITER))
        error("Newton's method did not converge.");

      // Store the result in sln_fine.
      sln_fine.copy(&Psi_prev_newton);

      // Visualize intermediate solutions and mesh during adaptivity.
      sprintf(title, "Solution, time level %d, adapt step %d", n, as);
      magview.set_title(title);
      AbsFilter mag(&Psi_prev_newton);
      magview.show(&mag);
      sprintf(title, "Fine mesh, time level %d, adapt step %d", n, as);
      ordview.set_title(title);
      ordview.show(rnls.get_space(0));

      // Calculate element errors and total error estimate.
      H1Adapt hp(&space);
      hp.set_solutions(&sln_coarse, &sln_fine);
      err_est = hp.calc_error() * 100;   // relative h1-error in percent
      info("ndof_coarse: %d, ndof_fine: %d, err_est: %g%%", 
	   space.get_num_dofs(), rnls.get_space(0)->get_num_dofs(), err_est);

      // If err_est too large, adapt the mesh.
      if (err_est < ERR_STOP) done = true;
      else {
        hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
        ndof = assign_dofs(&space);
        if (ndof >= NDOF_STOP) done = true;

        // Project the fine mesh solution on the new coarse mesh.
        info("---- Time step %d, adaptivity step %d, projecting fine mesh solution on new coarse mesh:",
          n, as);
        nls.project_global(&sln_fine, &Psi_prev_newton, PROJ_TYPE);

        if (NEWTON_ON_COARSE_MESH) {
          // Newton's loop on the coarse mesh.
          info("---- Time step %d, adaptivity step %d, Newton solve on new coarse mesh:", n, as);
          if (!nls.solve_newton(&Psi_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER))
            error("Newton's method did not converge.");
        }

        // Store the result in sln_coarse.
        sln_coarse.copy(&Psi_prev_newton);

        as++;
      }
    }
    while (!done);

    // Add entries to convergence graphs.
    graph_time_err.add_values(n*TAU, err_est);
    graph_time_err.save("time_error.dat");
    graph_time_dof.add_values(n*TAU, space.get_num_dofs());
    graph_time_dof.save("time_dof.dat");

    // Copy result of the Newton's iteration into Psi_prev_time.
    Psi_prev_time.copy(&sln_fine);
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
