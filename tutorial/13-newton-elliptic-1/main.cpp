#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function.h"

//  This example shows an introductory application of the Newton's
//  method to a nonlinear elliptic problem. We use zero Dirichlet boundary
//  conditions and a constant initial guess for the Newton's method.
//  The treatment of nonzero Dirichlet BC and a more general initial guess
//  will be shown in the next example newton-elliptic-2.
//
//  PDE: stationary heat transfer equation with nonlinear thermal
//  conductivity, - div[lambda(u)grad u] = 0
//
//  Domain: unit square (-10,10)^2
//
//  BC: Zero Dirichlet
//
//  The following parameters can be changed:

const int INIT_GLOB_REF_NUM = 3;      // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 5;       // Number of initial refinements towards boundary.
const int P_INIT = 2;                 // Initial polynomial degree.
const double NEWTON_TOL = 1e-6;       // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;      // Maximum allowed number of Newton iterations.
const double INIT_COND_CONST = 3.0;   // COnstant initial condition.

// Thermal conductivity (temperature-dependent)
// Note: for any u, this function has to be positive.
template<typename Real>
Real lam(Real u) { return 1 + pow(u, 4); }

// Derivative of the thermal conductivity with respect to 'u'.
template<typename Real>
Real dlam_du(Real u) { return 4*pow(u, 3); }

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int marker, double x, double y)
{
  return 0;
}

// Heat sources (can be a general function of 'x' and 'y').
template<typename Real>
Real heat_src(Real x, Real y)
{
  return 1.0;
}

// Initial condition. It will be projected on the FE mesh 
// to obtain initial coefficient vector for the Newton's method.
scalar init_cond(double x, double y, double& dx, double& dy)
{
  dx = 0;
  dy = 0;
  return INIT_COND_CONST;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1,INIT_BDY_REF_NUM);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Previous solution for the Newton's iteration.
  Solution u_prev;

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(jac), H2D_UNSYM, H2D_ANY, &u_prev);
  wf.add_vector_form(callback(res), H2D_ANY, &u_prev);

  // Initialize the linear system.
  NonlinSystem nls(&wf, &space);

  // Project the function init_cond() on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton'w method.");
  nls.project_global(init_cond, &u_prev);  

  // Show the initial condition for the Newton's method. 
  ScalarView pview("Projection of initial condition", 0, 0, 400, 300);
  pview.show(&u_prev);
  
  // Perform Newton's iteration.
  info("Performing Newton's iteration.");
  bool verbose = true; // Default is false.
  if (!nls.solve_newton(&u_prev, NEWTON_TOL, NEWTON_MAX_ITER, verbose)) 
    error("Newton's method did not converge.");

  // Visualise the solution and mesh.
  ScalarView sview("Solution", 410, 0, 400, 300);
  OrderView oview("Mesh", 820, 0, 400, 300);
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

