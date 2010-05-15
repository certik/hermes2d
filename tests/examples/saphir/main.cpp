#include "hermes2d.h"
#include "solver_umfpack.h"

using namespace RefinementSelectors;

//  This test makes sure that the example "saphir" works correctly.

const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;       // Number of initial uniform mesh refinements
const double THRESHOLD = 0.6;     // This is a quantitative parameter of the adapt(...) function and
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
const double ERR_STOP = 1e-2;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// Problem parameters
double LH = 96;                       // total horizontal length
double LH0 = 18;                      // first horizontal length
double LH1 = 48;                      // second horizontal length
double LH2 = 78;                      // third horizontal length
double LV = 96;                       // total vertical length
double LV0 = 18;                      // first vertical length
double LV1 = 48;                      // second vertical length
double LV2 = 78;                      // third vertical length
double Q_EXT = 1.0;                   // neutron source
double SIGMA_T_1 = 0.60;              // total cross-sections
double SIGMA_T_2 = 0.48;
double SIGMA_T_3 = 0.70;
double SIGMA_T_4 = 0.85;
double SIGMA_T_5 = 0.90;
double SIGMA_S_1 = 0.53;              // scattering cross sections
double SIGMA_S_2 = 0.20;
double SIGMA_S_3 = 0.66;
double SIGMA_S_4 = 0.50;
double SIGMA_S_5 = 0.89;
double Q_EXT_1 = 1;                   // sources
double Q_EXT_2 = 0;
double Q_EXT_3 = 1;
double Q_EXT_4 = 0;
double Q_EXT_5 = 0;

// Additional constants
double D_1 = 1/(3.*SIGMA_T_1);        //diffusion coefficients
double D_2 = 1/(3.*SIGMA_T_2);
double D_3 = 1/(3.*SIGMA_T_3);
double D_4 = 1/(3.*SIGMA_T_4);
double D_5 = 1/(3.*SIGMA_T_5);
double SIGMA_A_1 = SIGMA_T_1 - SIGMA_S_1;  // absorption coefficients
double SIGMA_A_2 = SIGMA_T_2 - SIGMA_S_2;
double SIGMA_A_3 = SIGMA_T_3 - SIGMA_S_3;
double SIGMA_A_4 = SIGMA_T_4 - SIGMA_S_4;
double SIGMA_A_5 = SIGMA_T_5 - SIGMA_S_5;

/********** Boundary conditions ***********/

// Boundary condition types
int bc_types(int marker)
{
  if (marker == 1) return BC_NATURAL;
  else return BC_ESSENTIAL;
}

// Dirichlet boundary condition values
scalar bc_values(int marker, double x, double y)
{
  return 0.0;
}

/********** Weak forms ***********/

// Bilinear form (material 1)
template<typename Real, typename Scalar>
Scalar bilinear_form_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return D_1 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v)
         + SIGMA_A_1 * int_u_v<Real, Scalar>(n, wt, u, v);
}

// Bilinear form (material 2)
template<typename Real, typename Scalar>
Scalar bilinear_form_2(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return D_2 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v)
         + SIGMA_A_2 * int_u_v<Real, Scalar>(n, wt, u, v);
}

// Bilinear form (material 3)
template<typename Real, typename Scalar>
Scalar bilinear_form_3(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return D_3 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v)
         + SIGMA_A_3 * int_u_v<Real, Scalar>(n, wt, u, v);
}

// Bilinear form (material 4)
template<typename Real, typename Scalar>
Scalar bilinear_form_4(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return D_4 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v)
         + SIGMA_A_4 * int_u_v<Real, Scalar>(n, wt, u, v);
}

// Bilinear form (material 5)
template<typename Real, typename Scalar>
Scalar bilinear_form_5(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return D_5 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v)
         + SIGMA_A_5 * int_u_v<Real, Scalar>(n, wt, u, v);
}

// Integration order for the bilinear forms
Ord bilinear_form_ord(int n, double *wt, Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return u->val[0] * v->val[0]; // returning the sum of the degrees of the basis
                                // and test function
}

// Linear form (material 1)
template<typename Real, typename Scalar>
Scalar linear_form_1(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_v<Real, Scalar>(n, wt, v);
}

// Linear form (material 3)
template<typename Real, typename Scalar>
Scalar linear_form_3(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_v<Real, Scalar>(n, wt, v);
}

// Integration order for the linear forms
Ord linear_form_ord(int n, double *wt, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return v->val[0];  // q_ext is piecewise constant, thus
                     // returning the polynomial degree of the test function;
}

int main(int argc, char* argv[])
{
  // Load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);
  // initial uniform mesh refinement
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // Create finite element space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // Enumerate basis functions
  space.assign_dofs();

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form_1, bilinear_form_ord, H2D_SYM, 1);
  wf.add_biform(0, 0, bilinear_form_2, bilinear_form_ord, H2D_SYM, 2);
  wf.add_biform(0, 0, bilinear_form_3, bilinear_form_ord, H2D_SYM, 3);
  wf.add_biform(0, 0, bilinear_form_4, bilinear_form_ord, H2D_SYM, 4);
  wf.add_biform(0, 0, bilinear_form_5, bilinear_form_ord, H2D_SYM, 5);
  wf.add_liform(0, linear_form_1, linear_form_ord, 1);
  wf.add_liform(0, linear_form_3, linear_form_ord, 3);

  // Matrix solver
  UmfpackSolver solver;

  // create a selector which will select optimal candidate
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER, &shapeset);

  // Convergence graph wrt. the number of degrees of freedom
  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Error Convergence for the Saphir Problem", "Degrees of Freedom", "Error Estimate [%]");
  graph.add_row("error estimate", "k", "-", "o");

  // Convergence graph wrt. CPU time
  GnuplotGraph graph_cpu;
  graph_cpu.set_captions("Error Convergence for the Saphir Problem", "CPU Time", "Error Estimate [%]");
  graph_cpu.add_row("error estimate", "k", "-", "o");
  graph_cpu.set_log_y();

  // Adaptivity loop
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

    // Solve the coarse mesh problem
    LinSystem ls(&wf, &solver);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln_coarse);

    // Solve the fine mesh problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &sln_fine);

    // Calculate error estimate wrt. fine mesh solution
    H1Adapt hp(&space);
    hp.set_solutions(&sln_coarse, &sln_fine);
    double err_est = hp.calc_error() * 100;

    // time measurement
    cpu_time.tick();

    // report results
    info("Estimate of error: %g%%", err_est);

    // Add entry to DOF convergence graph
    graph.add_values(0, space.get_num_dofs(), err_est);
    graph.save("conv_dof.gp");

    // Add entry to CPU convergence graph
    graph_cpu.add_values(0, cpu_time.accumulated(), err_est);
    graph_cpu.save("conv_cpu.gp");

    // If err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      ndof = assign_dofs(&space);
      if (ndof >= NDOF_STOP) done = true;
    }

    // Time measurement
    cpu_time.tick();
  }
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  int n_dof_allowed = 5850;
  printf("n_dof_actual = %d\n", ndof);
  printf("n_dof_allowed = %d\n", n_dof_allowed);// ndofs was 5701 at the time this test was created
  if (ndof <= n_dof_allowed) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }

}
