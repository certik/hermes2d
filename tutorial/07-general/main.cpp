#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

//  This example solves a general second-order linear equation with non-constant
//  coefficients, and shows how integration orders in linear and bilinear forms
//  can be defined manually.
//
//  PDE: -d/dx(a_11(x,y)du/dx) - d/dx(a_12(x,y)du/dy) - d/dy(a_21(x,y)du/dx) - d/dy(a_22(x,y)du/dy)
//       + a_1(x,y)du/dx + a_21(x,y)du/dy + a_0(x,y)u = rhs(x,y)
//
//  Domain: arbitrary
//
//  BC:  Dirichlet for boundary marker 1: u = g_D(x,y)
//       Natural for any other boundary marker:   (a_11(x,y)*nu_1 + a_21(x,y)*nu_2) * dudx
//                                              + (a_12(x,y)*nu_1 + s_22(x,y)*nu_2) * dudy = g_N(x,y)
//
//  The following parameters can be changed:

const int P_INIT = 2;             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 3;       // Number of initial uniform refinements.

// Problem parameters.
double a_11(double x, double y) {
  if (y > 0) return 1 + x*x + y*y;
  else return 1;
}

double a_22(double x, double y) {
  if (y > 0) return 1;
  else return 1 + x*x + y*y;
}

double a_12(double x, double y) {
  return 1;
}

double a_21(double x, double y) {
  return 1;
}

double a_1(double x, double y) {
  return 0.0;
}

double a_2(double x, double y) {
  return 0.0;
}

double a_0(double x, double y) {
  return 0.0;
}

double rhs(double x, double y) {
  return 1 + x*x + y*y;
}

double g_D(double x, double y) {
  return -cos(M_PI*x);
}

double g_N(double x, double y) {
  return 0;
}

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == 1) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return g_D(x, y);
}

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
  mloader.load("domain.mesh", &mesh);
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize the shapeset and the cache.
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // Create the finite element space.
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_essential_bc_values(essential_bc_values);
  space.set_uniform_order(P_INIT);

  // Enumerate degrees of freedom.
  int ndof = assign_dofs(&space);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_biform(bilinear_form, bilinear_form_ord, H2D_SYM);
  wf.add_liform(linear_form, linear_form_ord);
  wf.add_liform_surf(linear_form_surf, linear_form_surf_ord, 2);

  // Initialize views.
  ScalarView sview("Coarse solution", 0, 100, 798, 700);
  OrderView  oview("Polynomial orders", 800, 100, 798, 700);

  // Matrix solver.
  UmfpackSolver solver;

  // Solve the problem.
  Solution sln;
  LinSystem ls(&wf, &solver);
  ls.set_space(&space);
  ls.set_pss(&pss);
  ls.assemble();
  ls.solve(&sln);

  // Time measurement.
  cpu_time.tick();

  // View the solution and mesh.
  sview.show(&sln);
  oview.show(&space);

  // Skip visualization time.
  cpu_time.tick(H2D_SKIP);

  // Print timing information.
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
