#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#define DEBUG_ORDER
#include "hermes2d.h"
#include "solver_umfpack.h"

//  This example uses the Newton's method to solve a nonlinear complex-valued
//  time-dependent PDE (the Gross-Pitaevski equation describing the behavior
//  of Einstein-Bose quantum gases). For time-discretization one can use either
//  the first-order implicit Euler method or the second-order Crank-Nicolson
//  method.
//
//  PDE: non-stationary complex Gross-Pitaevski equation
//  describing resonances in Bose-Einstein condensates
//
//  ih \partial \psi/\partial t = -h^2/(2m) \Delta \psi +
//  g \psi |\psi|^2 + 1/2 m \omega^2 (x^2 + y^2) \psi
//
//  Domain: square (-1, 1)^2
//
//  BC:  homogeneous Dirichlet everywhere on the boundary

const int P_INIT = 4;            // Initial polynomial degree
const double TAU = 0.001;        // Time step
const double T_FINAL = 2;        // Time interval length
const int INIT_REF_NUM = 3;      // Number of initial uniform refinements
const int TIME_DISCR = 2;        // 1 for implicit Euler, 2 for Crank-Nicolson
const int PROJ_TYPE = 1;         // For the projection of the initial condition
                                 // on the initial mesh: 1 = H1 projection, 0 = L2 projection
const double NEWTON_TOL = 1e-5;  // Stopping criterion for the Newton's method
const int NEWTON_MAX_ITER = 100; // Maximum allowed number of Newton iterations

// Problem constants
const double H = 1;              // Planck constant 6.626068e-34;
const double M = 1;              // mass of boson
const double G = 1;              // coupling constant
const double OMEGA = 1;          // frequency

// Initial conditions
scalar fn_init(double x, double y, scalar& dx, scalar& dy)
{
  scalar val = exp(-10*(x*x + y*y));
  dx = val * (-20.0 * x);
  dy = val * (-20.0 * y);
  return val;
}

// Boundary condition types
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Boundary condition values
scalar bc_values(int marker, double x, double y)
{
 return 0;
}

// Weak forms
#include "forms.cpp"

// Implicit Euler method (1st-order in time)
template<typename Real, typename Scalar>
Scalar residuum_euler(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{  return F_euler(n, wt, v, e, ext);  }
template<typename Real, typename Scalar>
Scalar jacobian_euler(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{  return J_euler(n, wt, u, v, e, ext);  }

// Crank-Nicolson method (1st-order in time)
template<typename Real, typename Scalar>
Scalar residuum_cranic(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{  return F_cranic(n, wt, v, e, ext);  }
template<typename Real, typename Scalar>
Scalar jacobian_cranic(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{  return J_cranic(n, wt, u, v, e, ext);  }

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // initial mesh refinements
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

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

  // solutions for the Newton's iteration and time stepping
  Solution Psi_prev_time,
           Psi_prev_newton;

  // initialize the weak formulation
  WeakForm wf(1);
  if(TIME_DISCR == 1) {
    wf.add_biform(0, 0, callback(jacobian_euler), H2D_UNSYM, H2D_ANY, 1, &Psi_prev_newton);
    wf.add_liform(0, callback(residuum_euler), H2D_ANY, 2, &Psi_prev_newton, &Psi_prev_time);
  }
  else {
    wf.add_biform(0, 0, callback(jacobian_cranic), H2D_UNSYM, H2D_ANY, 1, &Psi_prev_newton);
    wf.add_liform(0, callback(residuum_cranic), H2D_ANY, 2, &Psi_prev_newton, &Psi_prev_time);
  }

  // initialize the nonlinear system and solver
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  // visualisation
  ScalarView view("", 0, 0, 700, 600);
  view.fix_scale_width(80);

  // setting initial condition at zero time level
  Psi_prev_time.set_exact(&mesh, fn_init);
  Psi_prev_newton.set_exact(&mesh, fn_init);
  nls.set_ic(&Psi_prev_newton, &Psi_prev_newton, PROJ_TYPE);

  // time stepping loop
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int n = 1; n <= nstep; n++)
  {

    info("---- Time step %d:", n);

    // Newton's method
    if (!nls.solve_newton_1(&Psi_prev_newton, NEWTON_TOL, NEWTON_MAX_ITER)) error("Newton's method did not converge.");

    // show the new time level solution
    char title[100];
    sprintf(title, "Time level %d", n);
    view.set_title(title);
    view.show(&Psi_prev_newton);

    // copy result of the Newton's iteration into Psi_prev_time
    Psi_prev_time.copy(&Psi_prev_newton);
  }

  // wait for all views to be closed
  View::wait();
  return 0;
}
