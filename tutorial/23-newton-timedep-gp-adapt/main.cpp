#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example shows how to combine automatic adaptivity with the Newton's
//  method for a nonlinear complex-valued time-dependent PDE (the Gross-Pitaevski
//  equation describing the behavior of Einstein-Bose quantum gases)
//  discretized implicitly in time (via implicit Euler or Crank-Nicolson).
//
//  PDE: non-stationary complex Gross-Pitaevski equation
//  describing resonances in Bose-Einstein condensates.
//
//  ih \partial \psi/\partial t = -h^2/(2m) \Delta \psi +
//  g \psi |\psi|^2 + 1/2 m \omega^2 (x^2 + y^2) \psi.
//
//  Domain: square (-1, 1)^2.
//
//  BC:  homogeneous Dirichlet everywhere on the boundary.
//
//  Time-stepping: either implicit Euler or Crank-Nicolson.
//
//  The following parameters can be changed:

const bool SOLVE_ON_COARSE_MESH = false;   // true... Newton is done on coarse mesh in every adaptivity step,
                                           // false...Newton is done on coarse mesh only once, then projection
                                           // of the fine mesh solution to coarse mesh is used.
const int INIT_REF_NUM = 2;                // Number of initial uniform refinements.
const int P_INIT = 2;                      // Initial polynomial degree.
const int TIME_DISCR = 2;                  // 1 for implicit Euler, 2 for Crank-Nicolson.
const double T_FINAL = 200.0;              // Time interval length.
const double TAU = 0.005;                  // Time step.

// Adaptivity.
const int UNREF_FREQ = 1;                  // Every UNREF_FREQ time step the mesh is unrefined.
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
const double ERR_STOP = 5.0;               // Stopping criterion for hp-adaptivity
                                           // (relative error between reference and coarse solution in percent)
const int NDOF_STOP = 60000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.

// Newton's method.
const double NEWTON_TOL_COARSE = 0.01;     // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.05;       // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 20;            // Maximum allowed number of Newton iterations.

// Problem parameters.
const double H = 1;                      // Planck constant 6.626068e-34.
const double M = 1;                      // Mass of boson.
const double G = 1;                      // Coupling constant.
const double OMEGA = 1;                  // Frequency.


// Initial condition.
scalar init_cond(double x, double y, scalar& dx, scalar& dy)
{
  scalar val = exp(-20*(x*x + y*y));
  dx = val * (-40.0*x);
  dy = val * (-40.0*y);
  return val;
}

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

// Weak forms.
# include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("square.mesh", &basemesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Solutions for the Newton's iteration and adaptivity.
  Solution Psi_prev_time, Psi_prev_newton;

  // Initialize the weak formulation.
  WeakForm wf;
  if(TIME_DISCR == 1) {
    wf.add_matrix_form(callback(jacobian_euler), H2D_UNSYM, H2D_ANY, &Psi_prev_newton);
    wf.add_vector_form(callback(residual_euler), H2D_ANY, Tuple<MeshFunction*>(&Psi_prev_newton, &Psi_prev_time));
  }
  else {
    wf.add_matrix_form(callback(jacobian_cranic), H2D_UNSYM, H2D_ANY, &Psi_prev_newton);
    wf.add_vector_form(callback(residual_cranic), H2D_ANY, Tuple<MeshFunction*>(&Psi_prev_newton, &Psi_prev_time));
  }

  // Initialize the nonlinear system.
  NonlinSystem nls(&wf, &space);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_time_dof, graph_time_err;

  // Project the function init_cond() on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton'w method.");
  Psi_prev_time.set_exact(&mesh, init_cond);             // Psi_prev_time set equal to init_cond().
  nls.project_global(&Psi_prev_time, &Psi_prev_newton);  // Initial vector calculated here.

  // Show the projection of the initial condition.
  char title[100];
  ScalarView magview("Projection of initial condition", 0, 0, 600, 500);
  magview.fix_scale_width(80);
  AbsFilter mag(&Psi_prev_newton);
  magview.show(&mag);
  OrderView ordview("Initial mesh", 610, 0, 600, 500);
  ordview.fix_scale_width(80);
  ordview.show(&space);
  
  // Newton's loop on the coarse mesh.
  info("Solving on coarse mesh.");
  bool verbose = true; // Default is false.
  if (!nls.solve_newton(&Psi_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
    error("Newton's method did not converge.");

  // Store the result in sln_coarse.
  Solution sln_coarse, sln_fine;
  sln_coarse.copy(&Psi_prev_newton);

  // Create selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, MAX_ORDER);

  // Time stepping loop.
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= nstep; ts++)
  {
    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("---- Time step %d, global derefinement.", ts);
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);

      // Project the fine mesh solution on the globally derefined mesh.
      info("---- Time step %d:", ts);
      if (SOLVE_ON_COARSE_MESH) 
        info("Projecting fine mesh solution to obtain initial vector on globally derefined mesh.");
      else 
        info("Projecting fine mesh solution on globally derefined mesh for error calculation.");
      nls.project_global(&sln_fine, &Psi_prev_newton);

      if (SOLVE_ON_COARSE_MESH) {
        // Newton's loop on the globally derefined mesh.
        info("Solving on globally derefined mesh.", ts);
        if (!nls.solve_newton(&Psi_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
          error("Newton's method did not converge.");
      }

      // Store the result in sln_coarse.
      sln_coarse.copy(&Psi_prev_newton);
    }

    // Adaptivity loop.
    bool done = false;
    double space_err_est;
    int as = 1;
    do
    {
      info("---- Time step %d, adaptivity step %d:", ts, as);

      // Initialize reference nonlinear system.
      RefSystem rnls(&nls);

      // Set initial condition for the Newton's method on the fine mesh.
      if (as == 1) {
        info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
        rnls.project_global(&sln_coarse, &Psi_prev_newton);
      }
      else {
        info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
        rnls.project_global(&sln_fine, &Psi_prev_newton);
      }

      // Newton's method on fine mesh.
      info("Solving on fine mesh.");
      if (!rnls.solve_newton(&Psi_prev_newton, NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose))
        error("Newton's method did not converge.");

      // Store the result in sln_fine.
      sln_fine.copy(&Psi_prev_newton);

      // Calculate element errors and total error estimate.
      info("Calculating error.");
      H1Adapt hp(&nls);
      hp.set_solutions(&sln_coarse, &sln_fine);
      space_err_est = hp.calc_error() * 100;   // relative h1-error in percent
      info("ndof_coarse: %d, ndof_fine: %d, space_err_est: %g%%", 
	   nls.get_num_dofs(), rnls.get_num_dofs(), space_err_est);

      // If space_err_est too large, adapt the mesh.
      if (space_err_est < ERR_STOP) done = true;
      else {
        info("Adapting coarse mesh.");
        done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (nls.get_num_dofs() >= NDOF_STOP) done = true;

        // Project the fine mesh solution on the new coarse mesh.
        if (SOLVE_ON_COARSE_MESH) 
          info("Projecting fine mesh solution to obtain initial vector on new coarse mesh.");
        else 
          info("Projecting fine mesh solution on coarse mesh for error calculation.");
        nls.project_global(&sln_fine, &Psi_prev_newton);

        if (SOLVE_ON_COARSE_MESH) {
          // Newton's loop on the coarse mesh.
          info("---- Time step %d, adaptivity step %d, solving on new coarse mesh.", ts, as);
          if (!nls.solve_newton(&Psi_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
            error("Newton's method did not converge.");
        }

        // Store the result in sln_coarse.
        sln_coarse.copy(&Psi_prev_newton);

        as++;
      }
    }
    while (!done);

    // Visualize solution and mesh.
    sprintf(title, "Solution, time level %d", ts);
    magview.set_title(title);
    AbsFilter mag(&sln_coarse);
    magview.show(&mag);
    sprintf(title, "Mesh, time level %d", ts);
    ordview.set_title(title);
    ordview.show(&space);

    // Add entries to convergence graphs.
    graph_time_err.add_values(ts*TAU, space_err_est);
    graph_time_err.save("time_error.dat");
    graph_time_dof.add_values(ts*TAU, nls.get_num_dofs());
    graph_time_dof.save("time_dof.dat");

    // Copy result of the Newton's iteration into Psi_prev_time.
    Psi_prev_time.copy(&sln_fine);
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
