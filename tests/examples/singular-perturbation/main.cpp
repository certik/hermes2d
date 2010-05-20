#include "hermes2d.h"
#include "solver_umfpack.h"

using namespace RefinementSelectors;

///  This test makes sure that the example "singular-perturbation" works correctly.
///
///  Parameters
///  - INIT_REF_NUM=1
///  - INIT_REF_NUM_BDY=3
///  - P_INIT=1
///  - THRESHOLD=0.3
///  - STRATEGY=0
///  - CAND_LIST=HP_ANISO
///  - MESH_REGULARITY=-1
///  - CONV_EXP=1.0
///  - ERR_STOP=0.001
///  - NDOF_STOP=100000
///
///  Results for given parameters
///  - DOFs: 3081

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
const double ERR_STOP = 0.001;   // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;     // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// problem constants
const double K_squared = 1e4;     // Equation parameter.

// boundary condition types
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// function values for Dirichlet boundary conditions
scalar bc_values(int marker, double x, double y)
{
  return 0;
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) + K_squared * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return K_squared * int_v<Real, Scalar>(n, wt, v);
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

  // enumerate basis functions
  space.assign_dofs();

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form), H2D_SYM);
  wf.add_liform(0, callback(linear_form));

  // matrix solver
  UmfpackSolver solver;

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof_est, graph_cpu_est;

  // create a selector which will select optimal candidate
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER, &shapeset);

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
    LinSystem ls(&wf, &solver);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln_coarse);

    // solve the fine mesh problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &sln_fine);

    // calculate error estimate wrt. fine mesh solution
    H1Adapt hp(&space);
    hp.set_solutions(&sln_coarse, &sln_fine);
    double err_est = hp.calc_error() * 100;

    // time measurement
    cpu_time.tick();

    // report results
    info("Estimate of error: %g%%", err_est);

    // add entries to DOF convergence graph
    graph_dof_est.add_values(space.get_num_dofs(), err_est);
    graph_dof_est.save("conv_dof_est.dat");

    // add entries to CPU convergence graph
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est);
    graph_cpu_est.save("conv_cpu_est.dat");

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
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  int n_dof_allowed = 3300;
  printf("n_dof_actual = %d\n", ndof);
  printf("n_dof_allowed = %d\n", n_dof_allowed);// ndofs was 625 at the time this test was created
  if (ndof <= n_dof_allowed) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
