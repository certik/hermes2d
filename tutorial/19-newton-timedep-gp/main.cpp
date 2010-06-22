#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#define DEBUG_ORDER
#include "hermes2d.h"

//  This example uses the Newton's method to solve a nonlinear complex-valued
//  time-dependent PDE (the Gross-Pitaevski equation describing the behavior
//  of Einstein-Bose quantum gases). For time-discretization one can use either
//  the first-order implicit Euler method or the second-order Crank-Nicolson
//  method.
//
//  PDE: non-stationary complex Gross-Pitaevski equation
//  describing resonances in Bose-Einstein condensates.
//
//  ih \partial \psi/\partial t = -h^2/(2m) \Delta \psi +
//  g \psi |\psi|^2 + 1/2 m \omega^2 (x^2 + y^2) \psi.
//
//  Domain: square (-1, 1)^2.
//
//  BC:  homogeneous Dirichlet everywhere on the boundary

const int INIT_REF_NUM = 2;      // Number of initial uniform refinements.
const int P_INIT = 4;            // Initial polynomial degree.
const double TAU = 0.005;        // Time step.
const double T_FINAL = 2;        // Time interval length.
const int TIME_DISCR = 2;        // 1 for implicit Euler, 2 for Crank-Nicolson.
const double NEWTON_TOL = 1e-5;  // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100; // Maximum allowed number of Newton iterations.

// Problem constants
const double H = 1;              // Planck constant 6.626068e-34.
const double M = 1;              // Mass of boson.
const double G = 1;              // Coupling constant.
const double OMEGA = 1;          // Frequency.

// Initial conditions.
scalar init_cond(double x, double y, scalar& dx, scalar& dy)
{
  scalar val = exp(-10*(x*x + y*y));
  dx = val * (-20.0 * x);
  dy = val * (-20.0 * y);
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
#include "forms.cpp"

// Implicit Euler method (1st-order in time).
template<typename Real, typename Scalar>
Scalar residual_euler(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{  return F_euler(n, wt, v, e, ext);  }
template<typename Real, typename Scalar>
Scalar jacobian_euler(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{  return J_euler(n, wt, u, v, e, ext);  }

// Crank-Nicolson method (1st-order in time).
template<typename Real, typename Scalar>
Scalar residual_cranic(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{  return F_cranic(n, wt, v, e, ext);  }
template<typename Real, typename Scalar>
Scalar jacobian_cranic(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{  return J_cranic(n, wt, u, v, e, ext);  }

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an H1 space.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Solutions for the Newton's iteration and time stepping.
  Solution Psi_prev_time,
           Psi_prev_newton;

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

  // Initialize views.
  ScalarView view("", 0, 0, 600, 500);
  view.fix_scale_width(80);

  // Project initial conditions on FE spaces to obtain initial coefficient 
  // vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton'w method.");
  Psi_prev_time.set_exact(&mesh, init_cond);              // Psi_prev_time is set equal to init_cond().
  nls.project_global(&Psi_prev_time, &Psi_prev_newton);   // Initial vector calculated here. 

  // Time stepping loop:
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= nstep; ts++)
  {

    info("---- Time step %d:", ts);

    // Newton's method.
    info("Performing Newton's method.");
    bool verbose = true; // Default is false.
    if (!nls.solve_newton(&Psi_prev_newton, NEWTON_TOL, NEWTON_MAX_ITER, verbose)) 
      error("Newton's method did not converge.");

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Time step %d", ts);
    view.set_title(title);
    view.show(&Psi_prev_newton);

    // Copy result of the Newton's iteration into Psi_prev_time.
    Psi_prev_time.copy(&Psi_prev_newton);
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
