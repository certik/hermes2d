#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

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

const int INIT_REF_NUM = 4;       // Number of initial uniform mesh refinements.
const int P_INIT = 3;             // Initial polynomial degree of all mesh elements.
const double NEWTON_TOL = 1e-6;   // Stopping criterion for the Newton's method
const int NEWTON_MAX_ITER = 100;  // Maximum allowed number of Newton iterations

const bool JFNK = false;          // true = jacobian-free method,
                                  // false = Newton
const int PRECOND = 2;            // Preconditioning by jacobian (1) or approximation of jacobian (2)
                                  // in case of jfnk,
                                  // default ML proconditioner in case of Newton

// Boundary condition types
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Exact solution and its derivatives
double exact(double x, double y, double &dx, double &dy)
{
	dx = (1- 2*x) * y * (1 - y);
	dy = (1- 2*y) * x * (1 - x);
	return  x * y * (1-x) * (1-y);
}

template<typename Real>
Real u(Real x, Real y)  {  return (x - x*x) * (y - y*y);  }
template<typename Real>
Real dudx(Real x, Real y)  {  return (1- 2*x) * y * (1 - y);  }
template<typename Real>
Real dudy(Real x, Real y)  {  return (1- 2*y) * x * (1 - x);  }

template<typename Real>
Real dudxx(Real x, Real y)  {  return -2.0 * (y-y*y);  }
template<typename Real>
Real dudyy(Real x, Real y)  {  return -2.0 * (x-x*x);  }
template<typename Real>
Real dudxy(Real x, Real y)  {  return (1- 2*y) * (1 - 2*x);  }

template<typename Real>
Real k(Real x, Real y)  {  return 1.0 / sqrt(1.0 + sqr(dudx(x,y)) + sqr(dudy(x,y)));  }
template<typename Real>
Real kx(Real x, Real y)  {  return -0.5 * pow(1.0 + sqr(dudx(x,y)) + sqr(dudy(x,y)), -1.5) *
                 (2.0 * dudx(x,y) * dudxx(x,y) + 2.0 * dudy(x,y) * dudxy(x,y));  }
template<typename Real>
Real ky(Real x, Real y)  {  return -0.5 * pow(1.0 + sqr(dudx(x,y)) + sqr(dudy(x,y)), -1.5) *
                 (2.0 * dudx(x,y) * dudxy(x,y) + 2.0 * dudy(x,y) * dudyy(x,y));  }

template<typename Real>
Real f(Real x, Real y)
{  return - kx(x,y) * dudx(x,y) - ky(x,y) * dudy(x,y) - k(x,y) * (dudxx(x,y) + dudyy(x,y)); }

// Weak forms
template<typename Real, typename Scalar>
Scalar jacobian_form_hermes(int n, double *wt, Func<Real> *vi, Func<Real> *vj, 
                            Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* u = ext->fn[0];
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( -0.5 * pow(1.0 + sqr(u->dx[i]) + sqr(u->dy[i]), -1.5) * 
                       (2.0 * u->dx[i] * vi->dx[i] + 2.0 * u->dx[i] * vi->dx[i])
                       * (u->dx[i] * vj->dx[i] + u->dy[i] * vj->dy[i]) +
                       (pow(1.0 + sqr(u->dx[i]) + sqr(u->dy[i]), -0.5))
                       * (vi->dx[i] * vj->dx[i] + vi->dy[i] * vj->dy[i]) );
  return result;
}

template<typename Real, typename Scalar>
Scalar residual_form_hermes(int n, double *wt, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* u = ext->fn[0];
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ((pow(1.0 + sqr(u->dx[i]) + sqr(u->dy[i]), -0.5)) * (u->dx[i] * vj->dx[i] + u->dy[i] * vj->dy[i])
                       - f(e->x[i], e->y[i]) * vj->val[i] );
  return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////
template<typename Real, typename Scalar>
Scalar jacobian_form_nox(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( -0.5 * pow(1.0 + sqr(u[0]->dx[i]) + sqr(u[0]->dy[i]), -1.5) *
                              (2.0 * u[0]->dx[i] * vi->dx[i] + 2.0 * u[0]->dx[i] * vi->dx[i])
                       * (u[0]->dx[i] * vj->dx[i] + u[0]->dy[i] * vj->dy[i]) +
                       (pow(1.0 + sqr(u[0]->dx[i]) + sqr(u[0]->dy[i]), -0.5))
                       * (vi->dx[i] * vj->dx[i] + vi->dy[i] * vj->dy[i]) );
  return result;
}

template<typename Real, typename Scalar>
Scalar precond_form_nox(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( vi->dx[i] * vj->dx[i] + vi->dy[i] * vj->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar residual_form_nox(int n, double *wt, Func<Real> *u[], Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ((pow(1.0 + sqr(u[0]->dx[i]) + sqr(u[0]->dy[i]), -0.5)) *
                        (u[0]->dx[i] * vj->dx[i] + u[0]->dy[i] * vj->dy[i])
                       - f(e->x[i], e->y[i]) * vj->val[i] );
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

  // Solutions
  Solution prev, sln1, sln2;

  // Initialize the shapeset and the cache,
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // Create an H1 space,
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_uniform_order(P_INIT);
  int ndof = assign_dofs(&space);
  info("Number of DOF: %d", space.get_num_dofs());

  // Time measurement,
  TimePeriod cpu_time;

  // First solve the problem using Hermes / NonlinSystem class / UMFpack matrix solver,

  info("!***************** Using NonlinSystem, Solving by Umfpack *****************");
  
  // Matrix solver,
  UmfpackSolver umfpack;

  // Define zero function on the mesh
  prev.set_zero(&mesh);

  // Initialize weak formulation,
  WeakForm wf1;
  wf1.add_biform(callback(jacobian_form_hermes), H2D_UNSYM, H2D_ANY, 1, &prev);
  wf1.add_liform(callback(residual_form_hermes), H2D_ANY, 1, &prev);

  // Initialize NonlinSystem,
  NonlinSystem nls(&wf1, &umfpack);
  nls.set_space(&space);
  nls.set_pss(&pss);

  // Project the function "prev" on the FE space "space",
  nls.project_global(&prev, &prev);

  // Perform Newton's iteration,
  info("NonlinSystem: Starting Newton's iteration... ");
  if (!nls.solve_newton(&prev, NEWTON_TOL, NEWTON_MAX_ITER)) 
    error("Newton's method did not converge.");
  info("Done.");

  // Storing the solution in "sln1"
  sln1.copy(&prev);

  // CPU time needed by UMFpack
  double umf_time = cpu_time.tick().last();

  info("!******************** Using FeProblem, Solving by NOX ************************");

  cpu_time.tick(H2D_SKIP);
 
  // Define zero function (again)
  prev.set_zero(&mesh);

  // Project the function "prev" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solution...");
  nls.project_global(&prev, &prev);
  info("Done.");

  // Get the coefficient vector.
  scalar *vec = nls.get_solution_vector();
  
  // Measure the projection time.
  double proj_time = cpu_time.tick().last();

  // Initialize the weak formulation for Trilinos.
  WeakForm wf2(1, JFNK ? true : false);
  if (!JFNK || (JFNK && PRECOND == 1)) wf2.add_jacform(callback(jacobian_form_nox), H2D_SYM);
  if (JFNK && PRECOND == 2) wf2.add_jacform(callback(precond_form_nox), H2D_SYM);
  wf2.add_resform(callback(residual_form_nox));

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
    vec = solver.get_solution();
    sln2.set_fe_solution(&space, &pss, vec);

    info("Number of nonlin iters: %d (norm of residual: %g)", solver.get_num_iters(), solver.get_residual());
    info("Total number of iters in linsolver: %d (achieved tolerance in the last step: %g)", 
         solver.get_num_lin_iters(), solver.get_achieved_tol());
  }
  else
    error("NOX failed.");

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
  view1.show(&sln1);
  ScalarView view2("Solution 2", 600, 0, 500, 400);
  view2.show(&sln2);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
