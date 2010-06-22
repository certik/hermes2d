#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function.h"

//  This example shows how to set an arbitrary initial guess for the
//  Newton's method, and nonzero Dirichlet boundary conditions.
//
//  PDE: stationary heat transfer equation with nonlinear thermal
//  conductivity, - div[lambda(u)grad u] = 0
//
//  Domain: unit square (-10,10)^2
//
//  BC: Dirichlet, see function dir_lift() below.
//
//  The following parameters can be changed:

const int P_INIT = 2;             // Initial polynomial degree
const double NEWTON_TOL = 1e-6;   // Stopping criterion for the Newton's method
const int NEWTON_MAX_ITER = 100;  // Maximum allowed number of Newton iterations
const int INIT_GLOB_REF_NUM = 3;  // Number of initial uniform mesh refinements
const int INIT_BDY_REF_NUM = 4;   // Number of initial refinements towards boundary

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
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1,INIT_BDY_REF_NUM);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Previous solution for the Newton's method.
  Solution u_prev;

  // Initialize the weak formulation
  WeakForm wf;
  wf.add_matrix_form(callback(jac), H2D_UNSYM, H2D_ANY, &u_prev);
  wf.add_vector_form(callback(res), H2D_ANY, &u_prev);

  // Initialize the linear system.
  NonlinSystem nls(&wf, &space);

  // Project the function init_cond() on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton'w method.");
  nls.project_global(init_cond, &u_prev);

  // Perform Newton's iteration.
  info("Performing Newton's iteration.");
  bool verbose = true; // Default is false.
  if (!nls.solve_newton(&u_prev, NEWTON_TOL, NEWTON_MAX_ITER, verbose)) 
    error("Newton's method did not converge.");

  // Visualise the solution and mesh.
  ScalarView sview("Solution", 0, 0, 400, 300);
  OrderView oview("Mesh", 410, 0, 400, 300);
  char title[100];
  sprintf(title, "Solution");
  sview.set_title(title);
  sview.show(&u_prev);
  sprintf(title, "Mesh");
  oview.set_title(title);
  oview.show(&space);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

