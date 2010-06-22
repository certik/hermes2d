#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function.h"

using namespace RefinementSelectors;

//  This example shows how the Newton's method can be combined with
//  automatic adaptivity.
//
//  PDE: stationary heat transfer equation with nonlinear thermal
//  conductivity, - div[lambda(u)grad u] = 0.
//
//  Domain: unit square (-10,10)^2.
//
//  BC: Dirichlet, see function dir_lift() below.
//
//  The following parameters can be changed:

const bool SOLVE_ON_COARSE_MESH = false;   // true...  Newton is done on coarse mesh in every adaptivity step.
                                           // false... Newton is done on coarse mesh only once, then projection
                                           // of the fine mesh solution to coarse mesh is used.
const int P_INIT = 1;                      // Initial polynomial degree.
const int INIT_GLOB_REF_NUM = 1;           // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 3;            // Number of initial refinements towards boundary.

const double THRESHOLD = 0.2;              // This is a quantitative parameter of the adapt(...) function and
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
const double ERR_STOP = 1.0;               // Stopping criterion for adaptivity (rel. error tolerance between the
                                           // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.
const double NEWTON_TOL_COARSE = 1e-4;     // Stopping criterion for the Newton's method on coarse mesh.
const double NEWTON_TOL_FINE = 1e-4;       // Stopping criterion for the Newton's method on fine mesh.
const int NEWTON_MAX_ITER = 100;           // Maximum allowed number of Newton iterations.

// Thermal conductivity (temperature-dependent).
// Note: for any u, this function has to be positive.
template<typename Real>
Real lam(Real u)
{
  return 1 + pow(u, 4);
}

// Derivative of the thermal conductivity with respect to 'u'.
template<typename Real>
Real dlam_du(Real u) {
  return 4*pow(u, 3);
}

// This function is used to define Dirichlet boundary conditions.
double dir_lift(double x, double y, double& dx, double& dy) {
  dx = (y+10)/10.;
  dy = (x+10)/10.;
  return (x+10)*(y+10)/100.;
}

// Initial condition. It will be projected on the FE mesh 
// to obtain initial coefficient vector for the Newton's method.
scalar init_cond(double x, double y, double& dx, double& dy)
{
  // Using the Dirichlet lift elevated by two
  double val = dir_lift(x, y, dx, dy) + 2;
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
  double dx, dy;
  return dir_lift(x, y, dx, dy);
}

// Heat sources (can be a general function of 'x' and 'y').
template<typename Real>
Real heat_src(Real x, Real y)
{
  return 1.0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Solutions for the Newton's iteration and adaptivity.
  Solution u_prev, sln_coarse, sln_fine;

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(jac), H2D_UNSYM, H2D_ANY, &u_prev);
  wf.add_vector_form(callback(res), H2D_ANY, &u_prev);

  // Initialize the coarse and fine mesh problems.
  NonlinSystem nls(&wf, &space);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu;

  // Project the function init_cond() on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector on coarse mesh.");
  nls.project_global(init_cond, &u_prev);

  // Initialize views.
  ScalarView sview_coarse("Coarse mesh solution", 0, 0, 350, 300); // coarse mesh solution
  OrderView oview_coarse("Coarse mesh", 360, 0, 350, 300);         // coarse mesh
  ScalarView sview_fine("Fine mesh solution", 720, 0, 350, 300);   // fine mesh solution
  OrderView oview_fine("Fine mesh", 1080, 0, 350, 300);            // fine mesh

  // Newton's loop on the coarse mesh.
  info("Solving on coarse mesh.");
  bool verbose = true; // Default is false.
  if (!nls.solve_newton(&u_prev, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) 
    error("Newton's method did not converge.");

  // Store the result in sln_coarse.
  sln_coarse.copy(&u_prev);

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Adaptivity loop:
  bool done = false; int as = 1;
  double err_est;
  do {
    info("---- Adaptivity step %d:", as);

    // Time measurement..
    cpu_time.tick();

    // Show coarse mesh and solution.
    sview_coarse.show(&sln_coarse);
    oview_coarse.show(&space);

    // Skip visualization.
    cpu_time.tick(H2D_SKIP);

    // Initialize the fine mesh problem.
    RefSystem rnls(&nls);

    // Set initial condition for the Newton's method on the fine mesh.
    if (as == 1) {
      info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
      rnls.project_global(&sln_coarse, &u_prev);
    }
    else {
      info("Projecting fine mesh solution to obtain initial vector on new fine mesh.");
      rnls.project_global(&sln_fine, &u_prev);
    }

    // Newton's loop on the fine mesh.
    info("Solving on fine mesh.");
    if (!rnls.solve_newton(&u_prev, NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose)) 
      error("Newton's method did not converge.");

    // Store the fine mesh solution in sln_fine.
    sln_fine.copy(&u_prev);

    // Time measurement.
    cpu_time.tick();

    // Show fine mesh and solution.
    sview_fine.show(&sln_fine);
    oview_fine.show(rnls.get_space(0));

    // Skip visualization time.
    cpu_time.tick(H2D_SKIP);

    // Calculate element errors and total error estimate.
    info("Calculating error.");
    H1Adapt hp(&nls);
    hp.set_solutions(&sln_coarse, &sln_fine);
    err_est = hp.calc_error() * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est: %g%%", 
      nls.get_num_dofs(), rnls.get_num_dofs(), err_est);

    // Add entry to DOF convergence graph.
    graph_dof.add_values(nls.get_num_dofs(), err_est);
    graph_dof.save("conv_dof.dat");

    // Add entry to CPU convergence graph.
    graph_cpu.add_values(cpu_time.accumulated(), err_est);
    graph_cpu.save("conv_cpu.dat");

    // If err_est too large, adapt the mesh.
    if (err_est < ERR_STOP) done = true;
    else {
      info("Adapting coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      if (nls.get_num_dofs() >= NDOF_STOP) {
        done = true;
        break;
      }

      // Project the fine mesh solution on the new coarse mesh.
      if (SOLVE_ON_COARSE_MESH) 
        info("Projecting fine mesh solution to obtain initial vector on new coarse mesh.");
      else 
        info("Projecting fine mesh solution on coarse mesh for error calculation.");
      nls.project_global(&sln_fine, &u_prev);

      if (SOLVE_ON_COARSE_MESH) {
        // Newton's loop on the new coarse mesh.
        info("Solving on coarse mesh.");
        if (!nls.solve_newton(&u_prev, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) 
          error("Newton's method did not converge.");
      }

      // Store the result in sln_coarse.
      sln_coarse.copy(&u_prev);
    }

    as++;
  }
  while (!done);

  // Time measurement.
  cpu_time.tick();

  verbose("Total running time: %g s", cpu_time.accumulated());

  // Wait for keyboard or mouse input.
  View::wait();
  return 0;
}

