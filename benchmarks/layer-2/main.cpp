#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

//  This singularly perturbed problem exhibits a thin boundary layer. The
//  exact solution facilitates convergence studies.
//
//  PDE: -Laplace u + K*K*u = K*K + g(x,y)
//
//  Domain: Square (-1,1)^2.
//
//  BC:  Homogeneous Dirichlet
//
//  Exact solution: v(x,y) = U(x)U(y) where U(t) = 1 - (exp(K*x)+exp(-K*x))/(exp(K) + exp(-K)) is
//  the exact solution to the 1D singularly perturbed problem -u'' + K*K*u = K*K* in (-1,1)
//  equipped with zero Dirichlet BC.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 1;       // Number of initial mesh refinements (the original mesh is just one element)
const int INIT_REF_NUM_BDY = 3;   // Number of initial mesh refinements towards the boundary
const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.3;     // This is a quantitative parameter of the adapt(...) function and
                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;           // Adaptive strategy:
                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                  //   error is processed. If more elements have similar errors, refine
                                  //   all to keep the mesh symmetric.
                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                  //   than THRESHOLD times maximum element error.
                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                  //   than THRESHOLD.
                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const RefinementSelectors::AllowedCandidates ADAPT_TYPE = RefinementSelectors::H2DRS_CAND_HP;         // Type of automatic adaptivity.
const bool ISO_ONLY = false;      // Isotropic refinement flag (concerns quadrilateral elements only).
                                  // ISO_ONLY = false ... anisotropic refinement of quad elements
                                  // is allowed (default),
                                  // ISO_ONLY = true ... only isotropic refinements of quad elements
                                  // are allowed.
const int MESH_REGULARITY = -1;   // Maximum allowed level of hanging nodes:
                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                  // Note that regular meshes are not supported, this is due to
                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;      // Default value is 1.0. This parameter influences the selection of
                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.0001;   // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;     // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// problem constants
const double K = 1e2;             // Equation parameter.

// solution to the 1D problem -u'' + K*K*u = K*K in (-1,1) with zero Dirichlet BC
double uhat(double x) {
  return 1. - (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K));
}
double duhat_dx(double x) {
  return -K * (exp(K*x) - exp(-K*x)) / (exp(K) + exp(-K));
}
double dduhat_dxx(double x) {
  return -K*K * (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K));
}

// exact solution u(x,y) to the 2D problem is defined as the
// Cartesian product of the 1D solutions
static double sol_exact(double x, double y, double& dx, double& dy)
{
  dx = duhat_dx(x) * uhat(y);
  dy = uhat(x) * duhat_dx(y);
  return uhat(x) * uhat(y);
}

// right-hand side
double rhs(double x, double y) {
  return -(dduhat_dxx(x)*uhat(y) + uhat(x)*dduhat_dxx(y)) + K*K*uhat(x)*uhat(y);
}

// boundary condition types
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// function values for Dirichlet boundary conditions
scalar bc_values(int marker, double x, double y)
{
  return 0;
}

// weak forms
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) + K*K * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real, Scalar>(n, wt, rhs, v, e);;
}

// integration order for the linear_form
Ord linear_form_ord(int n, double *wt, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(24);
}

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // initial mesh refinement
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, INIT_REF_NUM_BDY);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create finite element space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // enumerate degrees of freedom
  int ndof = assign_dofs(&space);

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form), H2D_SYM);
  wf.add_liform(0, linear_form, linear_form_ord);

  // visualize solution and mesh
  ScalarView sview("Coarse solution", 0, 0, 500, 400);
  OrderView  oview("Polynomial orders", 505, 0, 500, 400);

  // matrix solver
  UmfpackSolver solver;

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof_est, graph_dof_exact, graph_cpu_est, graph_cpu_exact;

  // selector
  RefinementSelectors::H1NonUniformHP ref_selector(ISO_ONLY, ADAPT_TYPE, CONV_EXP, H2DRS_DEFAULT_ORDER, &shapeset);

  // adaptivity loop
  int it = 1;
  bool done = false;
  TimePeriod cpu_time;
  Solution sln_coarse, sln_fine;
  do
  {
    info("!---- Adaptivity step %d ---------------------------------------------", it); it++;

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // solve the coarse mesh problem
    LinSystem ls(&wf, &solver);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln_coarse);

    // time measurement
    cpu_time.tick();

    // calculate error wrt. exact solution
    ExactSolution exact(&mesh, sol_exact);
    double error = h1_error(&sln_coarse, &exact) * 100;

    // view the solution
    sview.show(&sln_coarse);
    oview.show(&space);

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // solve the fine mesh problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &sln_fine);

    // calculate error estimate wrt. fine mesh solution
    H1AdaptHP hp(&space);
    double err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;

    // time measurement
    cpu_time.tick();

    // report results
    info("Exact solution error: %g%%", error);
    info("Estimate of error: %g%%", err_est);

    // add entries to DOF convergence graph
    graph_dof_exact.add_values(space.get_num_dofs(), error);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_dof_est.add_values(space.get_num_dofs(), err_est);
    graph_dof_est.save("conv_dof_est.dat");

    // add entries to CPU convergence graph
    graph_cpu_exact.add_values(cpu_time.accumulated(), error);
    graph_cpu_exact.save("conv_cpu_exact.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est);
    graph_cpu_est.save("conv_cpu_est.dat");

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(THRESHOLD, STRATEGY, &ref_selector, MESH_REGULARITY);
      ndof = assign_dofs(&space);
      if (ndof >= NDOF_STOP) done = true;
    }

    // time measurement
    cpu_time.tick();
  }
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // show the fine solution - this is the final result
  sview.set_title("Final solution");
  sview.show(&sln_fine);

  // wait for all views to be closed
  View::wait();
  return 0;
}
