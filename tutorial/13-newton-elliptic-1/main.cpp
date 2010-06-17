#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"
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

const int INIT_GLOB_REF_NUM = 3;  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 5;   // Number of initial refinements towards boundary.
const int P_INIT = 2;             // Initial polynomial degree.
const double NEWTON_TOL = 1e-6;   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;  // Maximum allowed number of Newton iterations.

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
  wf.add_biform(callback(jac), H2D_UNSYM, H2D_ANY, &u_prev);
  wf.add_liform(callback(res), H2D_ANY, &u_prev);

  // Matrix solver.
  UmfpackSolver solver;

  // Initialize the nonlinear system.
  NonlinSystem nls(&wf, &solver, &space);

  // Use a constant function as initial guess.
  double const_val = 3.0;
  u_prev.set_const(&mesh, const_val);

  // Project the function u_prev() on the FE space
  // to obtain initial guess u_prev for the Newton's method.
  info("Projecting initial condition on the FE space.");
  nls.project_global(&u_prev, &u_prev);

  // Perform Newton's iteration.
  info("Performing Newton's iteration.");
  if (!nls.solve_newton(&u_prev, NEWTON_TOL, NEWTON_MAX_ITER)) 
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

