#include "hermes2d.h"

using namespace RefinementSelectors;

// This test makes sure that example 10-adapt works correctly.

const int P_INIT = 2;                    // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.2;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
                                         // H2D_HP_ANISO_P, H2D_HP_ANISO. See User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double ERR_STOP = 1.0;             // Stopping criterion for adaptivity (rel. error tolerance between the
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows
                                         // over this limit. This is to prevent h-adaptivity to go on forever.

// Problem parameters.
const int OMEGA_1 = 1;
const int OMEGA_2 = 2;
const int STATOR_BDY = 2;
const double EPS_1 = 1.0;                // Relative electric permittivity in Omega_1.
const double EPS_2 = 10.0;               // Relative electric permittivity in Omega_2.
const double VOLTAGE = 50.0;             // Voltage on the stator.

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return (ess_bdy_marker == STATOR_BDY) ? VOLTAGE : 0.0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Check input parameters.
  // If true, coarse mesh FE problem is solved in every adaptivity step.
  // If false, projection of the fine mesh solution on the coarse mesh is used. 
  bool SOLVE_ON_COARSE_MESH = false;
  if (argc > 1 && strcmp(argv[1], "-coarse_mesh") == 0)
    SOLVE_ON_COARSE_MESH = true;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("motor.mesh", &mesh);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(biform1), H2D_SYM, OMEGA_1);
  wf.add_matrix_form(callback(biform2), H2D_SYM, OMEGA_2);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize the coarse mesh problem.
  LinSystem ls(&wf, &space);

  // Adaptivity loop:
  Solution sln_coarse, sln_fine;
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Assemble and solve the fine mesh problem.
    info("Solving on fine mesh.");
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(&sln_fine);    

    // Either solve on coarse mesh or project the fine mesh solution 
    // on the coarse mesh.
    if (SOLVE_ON_COARSE_MESH) {
      info("Solving on coarse mesh.");
      ls.assemble();
      ls.solve(&sln_coarse);
    }
    else {
      info("Projecting fine mesh solution on coarse mesh.");
      ls.project_global(&sln_fine, &sln_coarse);
    }

    // Calculate element errors and total error estimate.
    info("Calculating error.");
    H1Adapt hp(&ls);
    hp.set_solutions(&sln_coarse, &sln_fine);
    double err_est = hp.calc_error() * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est: %g%%", 
      ls.get_num_dofs(), rs.get_num_dofs(), err_est);

    // If err_est too large, adapt the mesh.
    if (err_est < ERR_STOP) done = true;
    else {
      info("Adapting coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      if (ls.get_num_dofs() >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());

  int ndof = ls.get_num_dofs();

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  if (ndof < 1100) {      // ndofs was 1025 at the time this test was created
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}

