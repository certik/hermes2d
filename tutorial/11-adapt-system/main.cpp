#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

using namespace RefinementSelectors;

// This example explains how to use the multimesh adaptive hp-FEM,
// where different physical fields (or solution components) can be
// approximated using different meshes and equipped with mutually
// independent adaptivity mechanisms. For the tutorial purposes,
// we manufactured an exact solution for a simplified version of
// the FitzHugh-Nagumo equation. This equation, in its full form,
// is a prominent example of activator-inhibitor systems in two-component
// reaction-diffusion equations, It describes a prototype of an
// excitable system (e.g., a neuron).
//
// PDE: Linearized FitzHugh-Nagumo equation
//      -d_u^2 \Delta u - f(u) + \sigma v = g_1,
//      -d_v^2 \Delta v - u + v = g_2.
// In the original equation, f(u) = \lambda u - u^3 - \kappa. For
// simplicity, here we just take f(u) = u.
//
// Domain: Square (-1,1)^2.
//
// BC: Both solution components are zero on the boundary.
//
// Exact solution: The functions g_1 and g_2 were calculated so that
//                 the exact solution is:
//        u(x,y) = U(x)*U(y) where U(t) = cos(M_PI*t/2)
//        v(x,y) = V(x)V(y) where V(t) = 1 - (exp(K*t)+exp(-K*t))/(exp(K) + exp(-K))
// Note: V(t) is the exact solution of the 1D singularly perturbed equation
//       -u'' + K*K*u = K*K in (-1, 1) with zero Dirichlet BC.
//
// The following parameters can be changed: In particular, compare hp- and
// h-adaptivity via the CAND_LIST option, and compare the multi-mesh vs.
// single-mesh using the MULTI parameter.

const int P_INIT_U = 2;          // Initial polynomial degree for u.
const int P_INIT_V = 2;          // Initial polynomial degree for v.
const int INIT_REF_BDY = 3;      // Number of initial boundary refinements
const bool MULTI = true;         // MULTI = true  ... use multi-mesh,
                                 // MULTI = false ... use single-mesh.
                                 // Note: In the single mesh option, the meshes are
                                 // forced to be geometrically the same but the
                                 // polynomial degrees can still vary.
const double THRESHOLD = 0.3;    // This is a quantitative parameter of the adapt(...) function and
                                 // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;          // Adaptive strategy:
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
                                         // See the Sphinx tutorial (http://hpfem.org/hermes2d/doc/src/tutorial-2.html#adaptive-h-fem-and-hp-fem) for details.
const int MESH_REGULARITY = -1;  // Maximum allowed level of hanging nodes:
                                 // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                 // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                 // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                 // Note that regular meshes are not supported, this is due to
                                 // their notoriously bad performance.
const double CONV_EXP = 1;       // Default value is 1.0. This parameter influences the selection of
                                 // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int MAX_ORDER = 10;        // Maximum allowed element degree
const double ERR_STOP = 0.5;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                 // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;     // Adaptivity process stops when the number of degrees of freedom grows over
                                 // this limit. This is mainly to prevent h-adaptivity to go on forever.

// problem constants
const double D_u = 1;
const double D_v = 1;
const double SIGMA = 1;
const double LAMBDA = 1;
const double KAPPA = 1;
const double K = 100;

// boundary condition types
BCType bc_types(int marker) { return BC_ESSENTIAL; }

// Dirichlet BC values
scalar bc_values(int marker, double x, double y) { return 0;}

// exact solution u(x,y) = U(x)*U(y) and its derivatives
double U(double t) {
  return cos(M_PI*t/2);
}
double dUdt(double t) {
  return -sin(M_PI*t/2)*(M_PI/2.);
}
double ddUdtt(double t) {
  return -cos(M_PI*t/2)*(M_PI/2.)*(M_PI/2.);
}
static double u_exact(double x, double y, double& dx, double& dy)
{
  dx = dUdt(x)*U(y);
  dy = U(x)*dUdt(y);
  return U(x)*U(y);
}

// exact solution v(x,y) = V(x)*V(y)
double V(double t) {
  return 1. - (exp(K*t) + exp(-K*t))/(exp(K) + exp(-K));
}
double dVdt(double t) {
  return -K*(exp(K*t) - exp(-K*t))/(exp(K) + exp(-K));
}
double ddVdtt(double t) {
  return -K*K*(exp(K*t) + exp(-K*t))/(exp(K) + exp(-K));
}
static double v_exact(double x, double y, double& dx, double& dy)
{
  dx = dVdt(x)*V(y);
  dy = V(x)*dVdt(y);
  return V(x)*V(y);
}

// right-hand side functions g_1 and g_2
double g_1(double x, double y)
{
  double Laplace_u = ddUdtt(x)*U(y) + U(x)*ddUdtt(y);
  double u = U(x)*U(y);
  double v = V(x)*V(y);
  return -D_u*D_u * Laplace_u - u + SIGMA*v;
}

double g_2(double x, double y)
{
  double Laplace_v = ddVdtt(x)*V(y) + V(x)*ddVdtt(y);
  double u = U(x)*U(y);
  double v = V(x)*V(y);
  return -D_v*D_v * Laplace_v - u + v;
}

// weak forms
#include "forms.cpp"

////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh umesh, vmesh;
  H2DReader mloader;
  mloader.load("square.mesh", &umesh);
  if (MULTI == false) umesh.refine_towards_boundary(1, INIT_REF_BDY);

  // create initial mesh (master mesh)
  vmesh.copy(&umesh);

  // initial mesh refinements in the vmesh towards the boundary
  if (MULTI == true) vmesh.refine_towards_boundary(1, INIT_REF_BDY);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset xpss(&shapeset);
  PrecalcShapeset ypss(&shapeset);

  // create the x displacement space
  H1Space uspace(&umesh, &shapeset);
  uspace.set_bc_types(bc_types);
  uspace.set_bc_values(bc_values);
  uspace.set_uniform_order(P_INIT_U);

  // create the y displacement space
  H1Space vspace(MULTI ? &vmesh : &umesh, &shapeset);
  vspace.set_bc_types(bc_types);
  vspace.set_bc_values(bc_values);
  vspace.set_uniform_order(P_INIT_V);

  // enumerate degrees of freedom
  int ndof = assign_dofs(2, &uspace, &vspace);

  // initialize the weak formulation
  WeakForm wf(2);
  wf.add_biform(0, 0, callback(bilinear_form_0_0));
  wf.add_biform(0, 1, callback(bilinear_form_0_1));
  wf.add_biform(1, 0, callback(bilinear_form_1_0));
  wf.add_biform(1, 1, callback(bilinear_form_1_1));
  wf.add_liform(0, linear_form_0, linear_form_0_ord);
  wf.add_liform(1, linear_form_1, linear_form_1_ord);

  // visualization of solution and meshes
  OrderView  uoview("Mesh for u", 0, 0, 500, 500);
  OrderView  voview("Mesh for v", 510, 0, 500, 500);
  ScalarView uview("Solution u", 1020, 0, 500, 500);
  ScalarView vview("Solution v", 1530, 0, 500, 500);

  // matrix solver
  UmfpackSolver umfpack;

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof, graph_cpu;

  // refinement selector
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER, &shapeset);

  // adaptivity loop
  int it = 1;
  bool done = false;
  TimePeriod cpu_time;
  Solution u_sln_coarse, v_sln_coarse;
  Solution u_sln_fine, v_sln_fine;
  do
  {
    info("---- Adaptivity step %d ---------------------------------------------", it); it++;

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // enumerate degrees of freedom
    ndof = assign_dofs(2, &uspace, &vspace);

    // solve the coarse mesh problem
    LinSystem ls(&wf, &umfpack);
    ls.set_spaces(2, &uspace, &vspace);
    ls.set_pss(2, &xpss, &ypss);
    ls.assemble();
    ls.solve(2, &u_sln_coarse, &v_sln_coarse);

    // view the solution and meshes
    cpu_time.tick();
    info("u_dof=%d, v_dof=%d", uspace.get_num_dofs(), vspace.get_num_dofs());
    uview.show(&u_sln_coarse);
    vview.show(&v_sln_coarse);
    uoview.show(&uspace);
    voview.show(&vspace);
    cpu_time.tick(H2D_SKIP);

    // solve the fine mesh problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(2, &u_sln_fine, &v_sln_fine);

    // calculate element errors and total error estimate
    H1Adapt hp(Tuple<Space*>(&uspace, &vspace));
    hp.set_solutions(Tuple<Solution*>(&u_sln_coarse, &v_sln_coarse), Tuple<Solution*>(&u_sln_fine, &v_sln_fine));
    hp.set_biform(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
    hp.set_biform(0, 1, bilinear_form_0_1<scalar, scalar>, bilinear_form_0_1<Ord, Ord>);
    hp.set_biform(1, 0, bilinear_form_1_0<scalar, scalar>, bilinear_form_1_0<Ord, Ord>);
    hp.set_biform(1, 1, bilinear_form_1_1<scalar, scalar>, bilinear_form_1_1<Ord, Ord>);
    double err_est = hp.calc_error(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_ABS) * 100;

    // time measurement
    cpu_time.tick();

    // calculate error wrt. exact solution
    ExactSolution uexact(&umesh, u_exact);
    ExactSolution vexact(&vmesh, v_exact);
    double u_error = h1_error(&u_sln_coarse, &uexact) * 100;
    double v_error = h1_error(&v_sln_coarse, &vexact) * 100;
    double error = std::max(u_error, v_error);

    // report results
    info("Exact solution error for u (H1 norm): %g%%", u_error);
    info("Exact solution error for v (H1 norm): %g%%", v_error);
    info("Exact solution error (maximum): %g%%", error);
    info("Estimate of error wrt. ref. solution (energy norm): %g%%", err_est);

    // add entry to DOF convergence graph
    graph_dof.add_values(uspace.get_num_dofs() + vspace.get_num_dofs(), error);
    if (MULTI == true) graph_dof.save("conv_dof_m.dat");
    else graph_dof.save("conv_dof_s.dat");

    // add entry to CPU convergence graph
    graph_cpu.add_values(cpu_time.accumulated(), error);
    if (MULTI == true) graph_cpu.save("conv_cpu_m.dat");
    else graph_cpu.save("conv_cpu_s.dat");

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // if err_est too large, adapt the mesh
    if (error < ERR_STOP) done = true;
    else {
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY, MULTI == true ? false : true);
      ndof = assign_dofs(2, &uspace, &vspace);
      if (ndof >= NDOF_STOP) done = true;
    }

    //time measurement
    cpu_time.tick();
  }
  while (!done);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // show the fine solution - this is the final result
  uview.set_title("Final solution u");
  uview.show(&u_sln_fine);
  vview.set_title("Final solution v");
  vview.show(&v_sln_fine);

  // wait for all views to be closed
  View::wait();
  return 0;
}
