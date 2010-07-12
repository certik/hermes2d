#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

//  The purpose of this example is to show how to use Trilinos for nonlinear PDE problems. It 
//  compares performance of The Newton's method in the NonlinSystem class in Hermes (using the 
//  Umfpack solver) with the performance of the Trilinos/NOX solver (using the Newton's method 
//  or JFNK, and with or without preconditioning).
//
//  PDE:  - \nabla (k \nabla u) = f
//  k = (1 + sqr(u_x) + sqr(u_y))^{-0.5}
//
//  Domain: Unit square.
//
//  BC: zero Dirichlet.
//
//  Exact solution: (x - x*x) * (y - y*y).
//
//  Initial guess for the Newton's method: zero function.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 4;       // Number of initial uniform mesh refinements.
const int P_INIT = 3;             // Initial polynomial degree of all mesh elements.
const double NEWTON_TOL = 1e-6;   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;  // Maximum allowed number of Newton iterations.

const bool JFNK = false;          // true = jacobian-free method,
                                  // false = Newton.
const int PRECOND = 2;            // Preconditioning by jacobian (1) or approximation of jacobian (2)
                                  // in case of JFNK,
                                  // Default ML proconditioner in case of Newton.

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Exact solution.
#include "exact_solution.cpp"

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
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Solutions.
  Solution prev, sln1, sln2;

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, NULL, P_INIT);
  info("Number of DOF: %d",  space.get_num_dofs());

  info("---- Using NonlinSystem, solving by Umfpack:");

  // Time measurement.
  cpu_time.tick(HERMES_SKIP);

  // Define zero function on the mesh
  prev.set_zero(&mesh);

  // Initialize weak formulation,
  WeakForm wf1;
  wf1.add_matrix_form(callback(jacobian_form_hermes), H2D_UNSYM, H2D_ANY, &prev);
  wf1.add_vector_form(callback(residual_form_hermes), H2D_ANY, &prev);

  // Initialize NonlinSystem,
  NonlinSystem nls(&wf1, &space);

  // Project the function "prev" on the FE space "space".
  info("Projecting initial condition on the FE space.");
  nls.project_global(&prev, &prev);

  // Perform Newton's iteration,
  info("Performing Newton's method.");
  if (!nls.solve_newton(&prev, NEWTON_TOL, NEWTON_MAX_ITER)) 
    error("Newton's method did not converge.");

  // Storing the solution in "sln1"
  sln1.copy(&prev);

  // CPU time needed by UMFpack
  double umf_time = cpu_time.tick().last();

  info("---- Using FeProblem, solving by NOX:");

  // Time measurement.
  cpu_time.tick(HERMES_SKIP);
 
  // Define zero function (again)
  prev.set_zero(&mesh);

  // Project the function "prev" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solution on the FE space.");
  nls.project_global(&prev, &prev);

  // Get the coefficient vector.
  scalar *vec = nls.get_solution_vector();
  
  // Measure the projection time.
  double proj_time = cpu_time.tick().last();

  // Initialize the weak formulation for Trilinos.
  WeakForm wf2(1, JFNK ? true : false);
  if (!JFNK || (JFNK && PRECOND == 1)) wf2.add_matrix_form(callback(jacobian_form_nox), H2D_SYM);
  if (JFNK && PRECOND == 2) wf2.add_matrix_form(callback(precond_form_nox), H2D_SYM);
  wf2.add_vector_form(callback(residual_form_nox));

  // Initialize FeProblem.
  H1Shapeset shapeset;
  FeProblem fep(&wf2);
  fep.set_spaces(&space);
  PrecalcShapeset pss(&shapeset);
  //fep.set_pss(1, &pss);

  // Initialize the NOX solver with the vector "vec".
  info("Initializing NOX.");
  NoxSolver solver(&fep);
  solver.set_init_sln(vec);

  // Choose preconditioning.
  MlPrecond pc("sa");
  if (PRECOND)
  {
    if (JFNK) solver.set_precond(&pc);
    else solver.set_precond("ML");
  }

  // Solve the matrix problem using NOX.
  info("Assembling by FeProblem, solving by NOX.");
  bool solved = solver.solve();
  if (solved)
  {
    vec = solver.get_solution_vector();
    sln2.set_fe_solution(&space, &pss, vec);

    info("Number of nonlin iterations: %d (norm of residual: %g)", solver.get_num_iters(), solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
         solver.get_num_lin_iters(), solver.get_achieved_tol());
  }
  else
    error("NOX failed.");

  // CPU time needed by NOX.
  double nox_time = cpu_time.tick().last();

  // Calculate exact errors.
  Solution ex;
  ex.set_exact(&mesh, &exact);
  info("Solution 1 (NonlinSystem - UMFpack): exact H1 error: %g (time %g s)", 
    100 * h1_error(&sln1, &ex), umf_time);
  info("Solution 2 (FeProblem - NOX):  exact H1 error: %g (time %g + %g s)", 
    100 * h1_error(&sln2, &ex), proj_time, nox_time);

  // Show both solutions.
  ScalarView view1("Solution 1", 0, 0, 500, 400);
  view1.show(&sln1);
  ScalarView view2("Solution 2", 600, 0, 500, 400);
  view2.show(&sln2);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
