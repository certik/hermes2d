#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

// This example shows how to run adaptive hp-FEM, h-FEM and p-FEM with
// basic control parameters. The underlying problem is a planar model
// of an electrostatic micromotor (MEMS). You may want to first run the model
// with hp-FEM (ADAPT_TYPE = 0), then with h-FEM (ADAPT_TYPE = 1), and
// last with p-FEM (ADAPT_TYPE = 2). The last option is there solely for
// completeness, since adaptive p-FEM is not really useful in practice.
// Uniform initial polynomial degree of mesh elements can be set using
// the variable P_INIT. The function adapt(...) takes the parameters THRESHOLD,
// STRATEGY, ADAPT_TYPE, ISO_ONLY, and MESH_REGULARITY whose meaning is
// explained below. Only the first parameter THRESHOLD is mandatory, all
// others can be left out and in that case default values are used.
// Additional control parameters are possible, these are demonstrated
// in the next tutorial example. Two types of convergence graphs are created
// -- error estimate wrt. the number of degrees of freedom (DOF), and error
// estimate wrt. CPU time. Later you will learn that also the error wrt.
// exact solution can be created for problems where exact solution is
// available. In this example you also can see how
// to define different material parameters in various parts of the
// computational domain.
//
// PDE: -div[eps_r(x,y) grad phi] = 0
//      eps_r = EPS1 in Omega_1 (surrounding air)
//      eps_r = EPS2 in Omega_2 (moving part of the motor)
//
// BC: phi = 0 V on Gamma_1 (left edge and also the rest of the outer boundary
//     phi = VOLTAGE on Gamma_2 (boundary of stator)
//
// The following parameters can be changed:

const int P_INIT = 2;             // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.2;     // This is a quantitative parameter of the adapt(...) function and
                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;           // Adaptive strategy:
                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                  //   error is processed. If more elements have similar errors, refine
                                  //   all to keep the mesh symmetric.
                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                  //   than THRESHOLD times maximum element error.
                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                  //   than THRESHOLD.
                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int ADAPT_TYPE = 0;         // Type of automatic adaptivity:
                                  // ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                  // ADAPT_TYPE = 1 ... adaptive h-FEM,
                                  // ADAPT_TYPE = 2 ... adaptive p-FEM.
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
const double ERR_STOP = 0.01;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// problem constants
const double EPS1 = 1.0;          // Relative electric permittivity in Omega_1.
const double EPS2 = 10.0;         // Relative electric permittivity in Omega_2.
const double VOLTAGE = 50.0;      // Voltage on the stator.

// boundary condition types
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Dirichlet boundary condition values
scalar bc_values(int marker, double x, double y)
{
  return (marker == 2) ? VOLTAGE : 0.0;
}

template<typename Real, typename Scalar>
Scalar biform1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return EPS1 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar biform2(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return EPS2 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("motor.mesh", &mesh);

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
  wf.add_biform(0, 0, callback(biform1), H2D_SYM, 1);
  wf.add_biform(0, 0, callback(biform2), H2D_SYM, 2);

  // visualize solution, gradient, and mesh
  ScalarView sview("Coarse solution", 0, 0, 600, 1000);
  VectorView gview("Gradient", 610, 0, 600, 1000);
  OrderView  oview("Polynomial orders", 1220, 0, 600, 1000);
  //gview.set_min_max_range(0.0, 400.0);

  // matrix solver
  UmfpackSolver solver;

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof, graph_cpu;

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

    // view the solution -- this can be slow; for illustration only
    sview.show(&sln_coarse);
    gview.show(&sln_coarse, &sln_coarse, H2D_EPS_NORMAL, H2D_FN_DX_0, H2D_FN_DY_0);
    oview.show(&space);

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // solve the fine mesh problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &sln_fine);

    // calculate element errors and total error estimate
    H1OrthoHP hp(1, &space);
    double err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;

    // time measurement
    cpu_time.tick();

    // report results
    info("Error estimate: %g%%", err_est);

    // add entry to DOF convergence graph
    graph_dof.add_values(space.get_num_dofs(), err_est);
    graph_dof.save("conv_dof.dat");

    // add entry to CPU convergence graph
    graph_cpu.add_values(cpu_time.accumulated(), err_est);
    graph_cpu.save("conv_cpu.dat");

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE, ISO_ONLY, MESH_REGULARITY);
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
  gview.show(&sln_fine, &sln_fine, H2D_EPS_HIGH, H2D_FN_DX_0, H2D_FN_DY_0);

  // wait for all views to be closed
  View::wait();
  return 0;
}
