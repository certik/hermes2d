#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

using namespace RefinementSelectors;

//  This example uses automatic adaptivity to solve a general second-order linear
//  equation with non-constant coefficients.
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

const bool SOLVE_ON_COARSE_MESH = false; // If true, coarse mesh FE problem is solved in every adaptivity step.
                                         // If false, projection of the fine mesh solution on the coarse mesh is used. 
const int P_INIT = 2;                    // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.6;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See the User Documentation for details.
const int MESH_REGULARITY = -1;   // Maximum allowed level of hanging nodes:
                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                  // Note that regular meshes are not supported, this is due to
                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;      // Default value is 1.0. This parameter influences the selection of
                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.01;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// boundary markers
const int BDY_DIRICHLET = 1;
const int BDY_NEUMANN = 2;

// Problem parameters
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

// Boundary condition types
BCType bc_types(int marker)
{
  if (marker == 1) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Dirichlet boundary condition values
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return g_D(x, y);
}

// (Volumetric) bilinear form
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    double x = e->x[i];
    double y = e->y[i];
    result += (a_11(x, y)*u->dx[i]*v->dx[i] +
               a_12(x, y)*u->dy[i]*v->dx[i] +
               a_21(x, y)*u->dx[i]*v->dy[i] +
               a_22(x, y)*u->dy[i]*v->dy[i] +
               a_1(x, y)*u->dx[i]*v->val[i] +
               a_2(x, y)*u->dy[i]*v->val[i] +
               a_0(x, y)*u->val[i]*v->val[i]) * wt[i];
  }
  return result;
}

// Integration order for the bilinear form
Ord bilinear_form_ord(int n, double *wt, Func<Ord> *u,
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return u->val[0] * v->val[0] * e->x[0] * e->x[0]; // returning the sum of the degrees of the basis
                                                    // and test function plus two
}

// Surface linear form (natural boundary conditions)
template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real, Scalar>(n, wt, g_N, v, e);
}

// Integration order for surface linear form
Ord linear_form_surf_ord(int n, double *wt, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return v->val[0] * e->x[0] * e->x[0];  // returning the polynomial degree of the test function plus two
}

// Volumetric linear form (right-hand side)
template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}

// Integration order for the volumetric linear form
Ord linear_form_ord(int n, double *wt, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return v->val[0] * e->x[0] * e->x[0];  // returning the polynomial degree of the test function plus two
}

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Refine the mesh uniformly.
  mesh.refine_all_elements();

  // Initialize the shapeset and the cache.
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // Create finite element space.
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
  wf.add_liform_surf(linear_form_surf, linear_form_surf_ord, BDY_NEUMANN);

  // Visualize solution and mesh.
  ScalarView sview("Coarse solution", 0, 0, 450, 350);
  OrderView  oview("Polynomial orders", 460, 0, 400, 350);

  // Matrix solver.
  UmfpackSolver solver;

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu;

  // Initialize refinement selector. 
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER, &shapeset);

  // Adaptivity loop:
  int as = 1;
  bool done = false;
  TimePeriod cpu_time;
  Solution sln_coarse, sln_fine;
  do
  {
    info("---- Adaptivity step %d:", as); as++;

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // Initialize the coarse mesh problem.
    LinSystem ls(&wf, &solver);
    ls.set_space(&space);
    ls.set_pss(&pss);

    // Initialize and solve the fine mesh problem.
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(&sln_fine);

    // Either solve on coarse mesh or project the fine mesh solution 
    // on the coarse mesh.
    if (SOLVE_ON_COARSE_MESH) {
      ls.assemble();
      ls.solve(&sln_coarse);
    }
    else ls.project_global(&sln_fine, &sln_coarse);

    // Time measurement.
    cpu_time.tick();

    // View the solution and mesh.
    sview.show(&sln_coarse);
    oview.show(&space);

    // Time measurement.
    cpu_time.tick(H2D_SKIP);

    // Calculate error estimate wrt. fine mesh solution.
    H1Adapt hp(&space);
    hp.set_solutions(&sln_coarse, &sln_fine);
    double err_est = hp.calc_error() * 100;

    // Time measurement.
    cpu_time.tick();

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est: %g%%", 
      space.get_num_dofs(), rs.get_space(0)->get_num_dofs(), err_est);

    // Add entry to DOF convergence graph.
    graph_dof.add_values(space.get_num_dofs(), err_est);
    graph_dof.save("conv_dof.dat");

    // Add entry to CPU convergence graph.
    graph_cpu.add_values(cpu_time.accumulated(), err_est);
    graph_cpu.save("conv_cpu.dat");

    // Time measurement.
    cpu_time.tick(H2D_SKIP);

    // If err_est too large, adapt the mesh.
    if (err_est < ERR_STOP) done = true;
    else {
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      ndof = assign_dofs(&space);
      if (ndof >= NDOF_STOP) done = true;
    }

    // Time measurement.
    cpu_time.tick();
  }
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the fine solution - the final result.
  sview.set_title("Final solution");
  sview.show(&sln_fine);

  // Wait for all views to be closed
  View::wait();
  return 0;
}
