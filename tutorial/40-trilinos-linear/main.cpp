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
const bool JFNK = true;          // true = jacobian-free method,
                                 // false = Newton
const bool PRECOND = true;       // Preconditioning by jacobian in case of jfnk,
                                 // default ML proconditioner in case of Newton

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Exact solution.
double exact(double x, double y, double &dx, double &dy)
{
	dx = 2*x;
	dy = 2*y;
	return x*x +y*y;
}

// Essential boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return x*x + y*y;
}

// Weak forms.
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ((-4.0) * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar jacobian_form(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( vi->dx[i] * vj->dx[i] + vi->dy[i] * vj->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar residual_form(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( u[0]->dx[i] * vi->dx[i] + u[0]->dy[i] * vi->dy[i] + 4.0 * vi->val[i] );
  return result;
}

int main(int argc, char **argv)
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();
 
  // Initialize shapeset and cache.
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // Create an H1 space.
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_essential_bc_values(essential_bc_values);
  space.set_uniform_order(P_INIT);

  // Enumerate degrees of freedom.
  int ndof = assign_dofs(&space);
  info("Number of DOF: %d", ndof);

  // Time measurement.
  TimePeriod cpu_time;

  // Solutions.
  Solution prev, sln1, sln2;

  info("!******************** Using LinSystem, Solving by UMFpack ************************");

  cpu_time.tick(H2D_SKIP);

  // Initialize matrix solver
  UmfpackSolver umfpack;

  // Initialize weak formulation.
  WeakForm wf1;
  wf1.add_biform(callback(bilinear_form));
  wf1.add_liform(callback(linear_form));

  // Initialize the linear system.
  LinSystem ls(&wf1, &umfpack);
  ls.set_space(&space);
  ls.set_pss(&pss);

  // Assemble and solve
  ls.assemble();
  ls.solve(&sln1);

  // CPU time needed by UMFpack
  double umf_time = cpu_time.tick().last();

  info("!******************** Using FeProblem, Solving by NOX Solver ************************");

  cpu_time.tick(H2D_SKIP);
 
  // Define zero function
  prev.set_zero(&mesh);
  
  // Project the function "prev" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solution...");
  ls.project_global(&prev, &prev);
  info("Done.");

  // Get the coefficient vector.
  scalar *vec = ls.get_solution_vector();
  
  // Measure the projection time.
  double proj_time = cpu_time.tick().last();
  
  // Initialize the weak formulation for Trilinos.
  WeakForm wf2(1, JFNK ? true : false);
  wf2.add_jacform(callback(jacobian_form), H2D_SYM);
  wf2.add_resform(callback(residual_form));

  // Initialize FeProblem.
  FeProblem fep(&wf2);
  fep.set_spaces(1, &space);
  fep.set_pss(1, &pss);

  // Initialize the NOX solver with the vector "vec".
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
  bool solved = solver.solve();
  if (solved)
  {
    double *s = solver.get_solution();
    sln2.set_fe_solution(&space, &pss, s);
    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      solver.get_num_iters(), solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      solver.get_num_lin_iters(), solver.get_achieved_tol());
  }
  else
    error("NOX failed");

  // CPU time needed by NOX
  double nox_time = cpu_time.tick().last();

  // Calculate exact errors
  Solution ex;
  ex.set_exact(&mesh, &exact);
  info("Solution 1 (Hermes-Newton): exact H1 error: %g (time %g s)", 
    100 * h1_error(&sln1, &ex), umf_time);
  info("Solution 2 (Trilinos-NOX):  exact H1 error: %g (time %g + %g s)", 
    100 * h1_error(&sln2, &ex), proj_time, nox_time);

  // Show both solutions
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
