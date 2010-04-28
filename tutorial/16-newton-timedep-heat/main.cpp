#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"
#include "function.h"

//  This example shows how the Newton's method is used for
//  a simple nonlinear parabolic PDE discretized in time
//  via the implicit Euler method.
//
//  PDE: time-dependent heat transfer equation with nonlinear thermal
//  conductivity, du/dt - div[lambda(u)grad u] = f
//
//  Domain: square (-10,10)^2
//
//  BC: Dirichlet, given by the function dir_lift() below.
//  IC: Same function dir_lift().
//
//  The following parameters can be changed:

const int P_INIT = 2;                  // Initial polynomial degree
const double TAU = 0.2;                // Time step
const double T_FINAL = 10.0;           // Time interval length
const int PROJ_TYPE = 1;               // For the projection of the initial condition
                                       // on the initial mesh: 1 = H1 projection, 0 = L2 projection
const int INIT_GLOB_REF_NUM = 3;       // Number of initial uniform mesh refinements
const int INIT_BDY_REF_NUM = 4;        // Number of initial refinements towards boundary

const double NEWTON_TOL = 1e-6;        // Stopping criterion for the Newton's method
const int NEWTON_MAX_ITER = 100;       // Maximum allowed number of Newton iterations

// Thermal conductivity (temperature-dependent)
// Note: for any u, this function has to be positive
template<typename Real>
Real lam(Real u)
{
  return 1 + pow(u, 4);
}

// Derivative of the thermal conductivity with respect to 'u'
template<typename Real>
Real dlam_du(Real u) {
  return 4*pow(u, 3);
}

// This function is used to define Dirichlet boundary conditions
double dir_lift(double x, double y, double& dx, double& dy) {
  dx = (y+10)/10.;
  dy = (x+10)/10.;
  return (x+10)*(y+10)/100.;
}

// Initial condition
scalar initial_condition(double x, double y, double& dx, double& dy)
{
  return dir_lift(x, y, dx, dy);
}

// Boundary condition type (essential = Dirichlet)
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Dirichlet boundary condition values
scalar bc_values(int marker, double x, double y)
{
  double dx, dy;
  return dir_lift(x, y, dx, dy);
}

// Heat sources (can be a general function of 'x' and 'y')
template<typename Real>
Real heat_src(Real x, Real y)
{
  return 1.0;
}

// Jacobian matrix
template<typename Real, typename Scalar>
Scalar jac(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev_newton = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] / TAU + dlam_du(u_prev_newton->val[i]) * u->val[i] *
                       (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
                       + lam(u_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
  return result;
}

// Fesidual vector
template<typename Real, typename Scalar>
Scalar res(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev_newton = ext->fn[0];
  Func<Scalar>* u_prev_time = ext->fn[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((u_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] / TAU +
                      lam(u_prev_newton->val[i]) * (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
		       - heat_src(e->x[i], e->y[i]) * v->val[i]);
  return result;
}

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // initial mesh refinements
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1,INIT_BDY_REF_NUM);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // enumerate degrees of freedom
  int ndof = assign_dofs(&space);

  // solutions for the Newton's iteration and
  // time stepping
  Solution u_prev_newton, u_prev_time;

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(jac), H2D_UNSYM, H2D_ANY, 1, &u_prev_newton);
  wf.add_liform(0, callback(res), H2D_ANY, 2, &u_prev_newton, &u_prev_time);

  // initialize the nonlinear system and solver
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  // project the function initial_condition() on the mesh
  nls.set_ic(initial_condition, &mesh, &u_prev_time, PROJ_TYPE);
  u_prev_newton.copy(&u_prev_time);

  // visualisation
  ScalarView sview("Solution", 0, 0, 500, 400);
  OrderView oview("Mesh", 520, 0, 450, 400);
  oview.show(&space);

  // time stepping loop
  double current_time = 0.0;
  int t_step = 1;
  do {
    info("---- Time step %d, t = %g s:", t_step, current_time); t_step++;

    // Newton's method
    if (!nls.solve_newton_1(&u_prev_newton, NEWTON_TOL, NEWTON_MAX_ITER)) error("Newton's method did not converge.");

    // update previous time level solution
    u_prev_time.copy(&u_prev_newton);

    // update time
    current_time += TAU;

    // show the new time level solution
    char title[100];
    sprintf(title, "Solution, t = %g", current_time);
    sview.set_title(title);
    sview.show(&u_prev_time);
  } while (current_time < T_FINAL);

  // wait for all views to be closed
  View::wait();
  return 0;
}

