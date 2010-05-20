#include "hermes2d.h"
#include "solver_umfpack.h"

using namespace RefinementSelectors;

/// This test makes sure that the benchmark "screen" works correctly.
///
///  Parameters
///  - P_INIT=1
///  - THERSHOLD=0.5
///  - STRATEGY=1
///  - CAND_LIST=HP_ANISO
///  - MESH_REGULARITY=-1
///  - ERR_STOP=0.1
///  - CONV_EXP=1.0
///  - NDOF_STOP=40000
///  - ERROR_WEIGHTS=(H: 1; P: 1; ANISO: 1)
///
///  Results for given parameters
///  - DOFs: 3994

const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.5;     // This is a quantitative parameter of the adapt(...) function and
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
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See the Sphinx tutorial (http://hpfem.org/hermes2d/doc/src/tutorial-2.html#adaptive-h-fem-and-hp-fem) for details.
const int MESH_REGULARITY = -1;   // Maximum allowed level of hanging nodes:
                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                  // Note that regular meshes are not supported, this is due to
                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;      // Default value is 1.0. This parameter influences the selection of
                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.1;      // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 40000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// problem constants
const double e_0  = 8.8541878176 * 1e-12;
const double mu_0 = 1.256 * 1e-6;
const double k = 1.0;

// exact solution
#include "exact_sol.cpp"

// boundary conditions
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

double2 tau[5] = { { 0, 0}, { 1, 0 },  { 0, 1 }, { -1, 0 }, { 0, -1 } };

cplx bc_values(int marker, double x, double y)
{
  scalar dx, dy;
  return exact0(x, y, dx, dy)*tau[marker][0] + exact1(x, y, dx, dy)*tau[marker][1];
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_curl_e_curl_f<Real, Scalar>(n, wt, u, v) - int_e_f<Real, Scalar>(n, wt, u, v);
}

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("screen-quad.mesh", &mesh);
//    mloader.load("screen-tri.mesh", &mesh);

  // initialize the shapeset and the cache
  HcurlShapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create finite element space
  HcurlSpace space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // enumerate basis functions
  space.assign_dofs();

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form), H2D_SYM);

  // matrix solver
  UmfpackSolver solver;

  // convergence graph wrt. the number of degrees of freedom
  GnuplotGraph graph;
  graph.set_captions("Error Convergence for the Screen Problem in H(curl)", "Degrees of Freedom", "Error [%]");
  graph.add_row("exact error", "k", "-", "o");
  graph.add_row("error estimate", "k", "--");
  graph.set_log_y();

  // convergence graph wrt. CPU time
  GnuplotGraph graph_cpu;
  graph_cpu.set_captions("Error Convergence for the Screen Problem in H(curl)", "CPU Time", "Error [%]");
  graph_cpu.add_row("exact error", "k", "-", "o");
  graph_cpu.add_row("error estimate", "k", "--");
  graph_cpu.set_log_y();

  // create a selector which will select optimal candidate
  HcurlProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER, &shapeset);
  selector.set_error_weights(1.0, 1.0, 1.0);

  // adaptivity loop
  int it = 1;
  int ndof;
  bool done = false;
  TimePeriod cpu_time;
  Solution sln_coarse, sln_fine;
  do
  {
    info("---- Adaptivity step %d ---------------------------------------------", it); it++;

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // solve the coarse mesh problem
    LinSystem sys(&wf, &solver);
    sys.set_spaces(1, &space);
    sys.set_pss(1, &pss);
    sys.assemble();
    sys.solve(1, &sln_coarse);

    // time measurement
    cpu_time.tick();

    // calculating error wrt. exact solution
    Solution ex;
    ex.set_exact(&mesh, exact);
    double error = 100 * hcurl_error(&sln_coarse, &ex);
    info("Exact solution error: %g%%", error);

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // solve the fine mesh problem
    RefSystem ref(&sys);
    ref.assemble();
    ref.solve(1, &sln_fine);

    // calculate error estimate wrt. fine mesh solution
    HcurlAdapt hp(&space);
    hp.set_solutions(&sln_coarse, &sln_fine);
    double err_est = hp.calc_error() * 100;

    // time measurement
    cpu_time.tick();

    // report results
    info("Error estimate: %g%%", err_est);

    // add entry to DOF convergence graph
    graph.add_values(0, space.get_num_dofs(), error);
    graph.add_values(1, space.get_num_dofs(), err_est);
    graph.save("conv_dof.gp");

    // add entry to CPU convergence graph
    graph_cpu.add_values(0, cpu_time.accumulated(),error);
    graph_cpu.add_values(1, cpu_time.accumulated(),err_est);
    graph_cpu.save("conv_cpu.gp");

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      ndof = assign_dofs(&space);
      if (ndof >= NDOF_STOP) done = true;
    }

    // time measurement
    cpu_time.tick();
  }
  while (!done);
  verbose("Total running time: %g s", cpu_time.accumulated());

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  int n_dof_allowed = 4000;
  printf("n_dof_actual = %d\n", ndof);
  printf("n_dof_allowed = %d\n", n_dof_allowed);// ndofs was 3161 at the time this test was created
  if (ndof <= n_dof_allowed) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}

