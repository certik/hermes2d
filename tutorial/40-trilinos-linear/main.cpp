#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

//  The purpose of this example is to show how to use Trilinos
//  for linear PDE problems. It compares performance of the LinSystem 
//  class in Hermes using the UMFpack matrix solver with the performance
//  of the Trilinos NOX solver (using Newton's method or JFNK, with or 
//  without preconditioning).
//
//  PDE: Poisson equation.
//
//  Domain: Square (-1, 1)^2.
//
//  BC: Nonhomogeneous Dirichlet, see the function essential_bc_values() below.
//
//  Exact solution: sqr(x) + sqr(y).
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 4;      // Number of initial uniform mesh refinements.
const int P_INIT = 2;            // Initial polynomial degree of all mesh elements.
const bool JFNK = true;          // true = Jacobian-free method,
                                 // false = Newton.
const bool PRECOND = true;       // Preconditioning by jacobian in case of jfnk,
                                 // default ML preconditioner in case of Newton.

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return x*x + y*y;
}

// Exact solution.
double exact(double x, double y, double &dx, double &dy)
{
	dx = 2*x;
	dy = 2*y;
	return x*x +y*y;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char **argv)
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();
 
  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  info("Number of DOF: %d", space.get_num_dofs());

  // Solutions.
  Solution prev, sln1, sln2;

  info("---- Using LinSystem, solving by UMFpack:");

  // Time measurement.
  cpu_time.tick(H2D_SKIP);

  // Matrix solver.
  UmfpackSolver umfpack;

  // Initialize weak formulation.
  WeakForm wf1;
  wf1.add_matrix_form(callback(bilinear_form));
  wf1.add_vector_form(callback(linear_form));

  // Initialize the linear system.
  LinSystem ls(&wf1, &umfpack, &space);

  // Assemble and solve.
  info("Assembling by LinSystem, solving by UMFpack.");
  ls.assemble();
  ls.solve(&sln1);

  // CPU time needed by UMFpack.
  double umf_time = cpu_time.tick().last();

  info("---- Using FeProblem, solving by NOX:");

  // Time measurement.
  cpu_time.tick(H2D_SKIP);
 
  // Define zero function.
  prev.set_zero(&mesh);
  
  // Project the function "prev" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solution on the FE mesh.");
  ls.project_global(&prev, &prev);

  // Get the coefficient vector.
  scalar *vec = ls.get_solution_vector();
  
  // Measure the projection time.
  double proj_time = cpu_time.tick().last();
  
  // Initialize the weak formulation for Trilinos.
  WeakForm wf2(1, JFNK ? true : false);
  wf2.add_jacform(callback(jacobian_form), H2D_SYM);
  wf2.add_resform(callback(residual_form));

  // FIXME: The entire FeProblem should be removed
  // and functionality merged into LinSystem.
  // Initialize FeProblem.
  H1Shapeset shapeset;
  FeProblem fep(&wf2);
  fep.set_spaces(&space);
  PrecalcShapeset pss(&shapeset);
  //fep.set_pss(Tuple<PrecalcShapeset*>(&pss));

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
    double *s = solver.get_solution_vector();
    sln2.set_fe_solution(&space, &pss, s);
    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      solver.get_num_iters(), solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      solver.get_num_lin_iters(), solver.get_achieved_tol());
  }
  else
    error("NOX failed");

  // CPU time needed by NOX.
  double nox_time = cpu_time.tick().last();

  // Calculate exact errors.
  Solution ex;
  ex.set_exact(&mesh, &exact);
  info("Solution 1 (LinSystem - UMFpack): exact H1 error: %g (time %g s)", 
    100 * h1_error(&sln1, &ex), umf_time);
  info("Solution 2 (FeProblem - NOX):  exact H1 error: %g (time %g + %g s)", 
    100 * h1_error(&sln2, &ex), proj_time, nox_time);

  // Show both solutions.
  ScalarView view1("Solution 1", 0, 0, 500, 400);
  view1.set_min_max_range(0, 2);
  view1.show(&sln1);
  ScalarView view2("Solution 2", 600, 0, 500, 400);
  view2.set_min_max_range(0, 2);
  view2.show(&sln2);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
