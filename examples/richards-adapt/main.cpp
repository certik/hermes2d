#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example uses adaptivity with dynamical meshes to solve
//  a simple version of the time-dependent Richard's equation.
//  The time discretization is backward Euler, and the Newton's
//  method is applied to solve the nonlinear problem in each time 
//  step. 
//
//  PDE: C(h)dh/dt - div(K(h)grad(h)) - (dK/dh)*(dh/dy) = 0
//  where K(h) = K_S*exp(alpha*h)                          for h < 0,
//        K(h) = K_S                                       for h >= 0,
//        C(h) = alpha*(theta_s - theta_r)*exp(alpha*h)    for h < 0,
//        C(h) = alpha*(theta_s - theta_r)                 for h >= 0.
//
//  Known exact solution, see the file exact_solution.cpp.
//
//  Domain: square (0, 100)^2.
//
//  BC: Dirichlet, given by the function init_cond() below.
//  IC: See init_cond().
//
//  The following parameters can be changed:

const double TIME_INIT = 5e-3;             // Initial time.
const int INIT_REF_NUM = 0;                // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 2;            // Number of initial mesh refinements towards the top edge.
const int P_INIT = 2;                      // Initial polynomial degree of all mesh elements.
const double TAU = 0.001;                  // Time step.
const double T_FINAL = 5.0;                // Time interval length.

// Adaptivity
const int UNREF_FREQ = 1;                  // Every UNREF_FREQth time step the mesh is unrefined.
const double THRESHOLD = 0.3;              // This is a quantitative parameter of the adapt(...) function and
                                           // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                    // Adaptive strategy:
                                           // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                           //   error is processed. If more elements have similar errors, refine
                                           //   all to keep the mesh symmetric.
                                           // STRATEGY = 1 ... refine all elements whose error is larger
                                           //   than THRESHOLD times maximum element error.
                                           // STRATEGY = 2 ... refine all elements whose error is larger
                                           //   than THRESHOLD.
                                           // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO_H;   // Predefined list of element refinement candidates. Possible values are
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
const double ERR_STOP = 0.1;               // Stopping criterion for adaptivity (rel. error tolerance between the
                                           // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.

// Newton's method
const double NEWTON_TOL_COARSE = 0.01;     // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.05;       // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 20;            // Maximum allowed number of Newton iterations.


// For the definition of initial condition.
int Y_POWER = 50;

// Problem parameters.
double K_S = 20.464;
double ALPHA = 1e-3;
double THETA_R = 0;
double THETA_S = 0.45;
double H_R = -1000;
double A = 100;
double L = 100;

// Current time.
double TIME = TIME_INIT;

// K:
double K(double h)
{
  if (h < 0) return K_S*exp(ALPHA*h);
  else return K_S;    
}

// dK/dh:
double dKdh(double h)
{
  if (h < 0) return K_S*ALPHA*exp(ALPHA*h);
  else return 0;
}

// ddK/dhh:
double ddKdhh(double h)
{
  if (h < 0) return K_S*ALPHA*ALPHA*exp(ALPHA*h);
  else return 0;
}

// C:
double C(double h)
{
  if (h < 0) return ALPHA*(THETA_S - THETA_R)*exp(ALPHA*h);
  else return ALPHA*(THETA_S - THETA_R);    
}

// dC/dh:
double dCdh(double h)
{
  if (h < 0) return ALPHA*(THETA_S - THETA_R)*ALPHA*exp(ALPHA*h);
  else return 0;    
}

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Exact solution.
#include "exact_solution.cpp"

// Essential (Dirichlet) boundary condition markers.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  double dx, dy;
  return exact_sol(x, y, dx, dy);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  /*
  for (int i=0; i<101; i++) {
    double y = i;
    double x = 50;
    double dx, dy;
    double h = exact_sol(x, y, dx, dy);
    printf("x = %g, y = %g, h = %g\n", x, y, h);
  }
  //exit(0);
  */

  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("square.mesh", &basemesh);

  // Perform initial mesh refinements.
  mesh.copy(&basemesh);
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(2, INIT_REF_NUM_BDY);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Adapt mesh to represent initial condition with given accuracy.
  bool verbose = true; 
  adapt_to_exact_function_h1(&space, exact_sol, &selector, THRESHOLD, STRATEGY, 
                             MESH_REGULARITY, ERR_STOP, NDOF_STOP, verbose);   
  info("Initial mesh: ndof = %d", space.get_num_dofs());

  // Solutions for the time stepping and the Newton's method.
  Solution u_prev_time, u_prev_newton;

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(jac, jac_ord, H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&u_prev_newton, &u_prev_time));
  wf.add_vector_form(res, res_ord, H2D_ANY, Tuple<MeshFunction*>(&u_prev_newton, &u_prev_time));
  
  // Initialize the nonlinear system.
  NonlinSystem nls(&wf, &space);

  // Error estimate and discrete problem size as a function of physical time.
  SimpleGraph graph_time_err_est, graph_time_dof_est, graph_time_err_exact;

  // Project the function init_cond() on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton'w method.");
  u_prev_time.set_exact(&mesh, exact_sol);            // u_prev_time set equal to exact_sol().
  nls.project_global(&u_prev_time, &u_prev_newton);   // Initial vector calculated here.

  // View the projection of the initial condition.
  ScalarView view("Projection of initial condition", 0, 0, 410, 300);
  OrderView ordview("Initial mesh", 420, 0, 350, 300);
  view.fix_scale_width(80);
  view.show(&u_prev_newton);
  ordview.show(&space);
  View::wait(H2DV_WAIT_KEYPRESS);

  // Newton's loop on the coarse mesh.
  info("Solving on coarse mesh.");
  verbose = true; // Default is false.
  if (!nls.solve_newton(&u_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
    error("Newton's method did not converge.");

  // Store the result in sln_coarse.
  Solution sln_coarse, sln_fine;
  sln_coarse.copy(&u_prev_newton);

  // Time stepping loop.
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    // Updating current time.
    TIME = ts*TAU;

    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);

      // Project fine mesh solution on the globally derefined mesh.
      info("---- Time step %d:", ts);
      info("Projecting fine mesh solution on globally derefined mesh for error calculation.");
      nls.project_global(&sln_fine, &u_prev_newton);

      // Store the result in sln_coarse.
      sln_coarse.copy(&u_prev_newton);
    }

    // Adaptivity loop (in space):
    bool done = false;
    double space_err_est, space_err_exact;
    int as = 1;
    do
    {
      info("---- Time step %d, adaptivity step %d:", ts, as);

      // Initialize reference nonlinear system.
      RefSystem rnls(&nls);

      // Set initial condition for the Newton's method on the fine mesh.
      if (as == 1) {
        info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
        rnls.project_global(&sln_coarse, &u_prev_newton);
      }
      else {
        info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
        rnls.project_global(&sln_fine, &u_prev_newton);
      }

      // Newton's method on fine mesh
      info("Solving on fine mesh.");
      if (!rnls.solve_newton(&u_prev_newton, NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose))
        error("Newton's method did not converge.");

      // Store the result in sln_fine.
      sln_fine.copy(&u_prev_newton);

      // Calculate error estimate wrt. fine mesh solution.
      info("Calculating error (est).");
      H1Adapt hp(&nls);
      hp.set_solutions(&sln_coarse, &sln_fine);
      space_err_est = hp.calc_error() * 100;

      // Calculate error wrt. exact solution.
      info("Calculating error (exact).");
      ExactSolution exact(&mesh, exact_sol);
      space_err_exact = h1_error(&sln_coarse, &exact) * 100;    

      info("ndof_coarse: %d, ndof_fine: %d, space_err_est: %g%%, space_err_exact: %g%%", 
	   nls.get_num_dofs(), rnls.get_num_dofs(), space_err_est, space_err_exact);

      // If space_err_est too large, adapt the mesh.
      if (space_err_est < ERR_STOP) done = true;
      else {
        info("Adapting coarse mesh.");
        done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (nls.get_num_dofs() >= NDOF_STOP) {
          done = true;
          break;
        }

        // Project the fine mesh solution on the new coarse mesh.
        info("Projecting fine mesh solution on coarse mesh for error calculation.");
        nls.project_global(&sln_fine, &u_prev_newton);

        // Store the result in sln_coarse.
        sln_coarse.copy(&u_prev_newton);

        as++;
      }
    }
    while (!done);

    // Visualize the solution and mesh.
    char title[100];
    sprintf(title, "Solution, time level %d", ts);
    view.set_title(title);
    view.show(&sln_coarse);
    sprintf(title, "Mesh, time level %d", ts);
    ordview.set_title(title);
    ordview.show(&space);

    // Add entries to convergence graphs.
    graph_time_err_est.add_values(ts*TAU, space_err_est);
    graph_time_err_est.save("time_error_est.dat");
    graph_time_dof_est.add_values(ts*TAU, nls.get_num_dofs());
    graph_time_dof_est.save("time_dof_est.dat");
    graph_time_err_exact.add_values(ts*TAU, space_err_exact);
    graph_time_err_exact.save("time_error_exact.dat");

    // Copy new time level solution into u_prev_time.
    u_prev_time.copy(&sln_fine);
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
