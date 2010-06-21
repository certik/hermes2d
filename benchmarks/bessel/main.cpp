#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example comes with an exact solution, and it describes the diffraction
//  of an electromagnetic wave from a re-entrant corner. Convergence graphs saved
//  (both exact error and error estimate, and both wrt. dof number and cpu time).
//
//  PDE: time-harmonic Maxwell's equations
//
//  Known exact solution, see functions exact_sol_val(), exact_sol(), exact()
//
//  Domain: L-shape domain
//
//  Meshes: you can use either "lshape3q.mesh" (quadrilateral mesh) or
//          "lshape3t.mesh" (triangular mesh). See the mesh.load(...) command below.
//
//  BC: perfect conductor on boundary markers 1 and 6 (essential BC)
//      impedance boundary condition on rest of boundary (natural BC)
//
//  The following parameters can be changed:

const bool SOLVE_ON_COARSE_MESH = false; // If true, coarse mesh FE problem is solved in every adaptivity step.
                                         // If false, projection of the fine mesh solution on the coarse mesh is used. 
const int INIT_REF_NUM = 1;              // Number of initial uniform mesh refinements.
const int P_INIT = 2;                    // Initial polynomial degree. NOTE: The meaning is different from
                                         // standard continuous elements in the space H1. Here, P_INIT refers
                                         // to the maximum poly order of the tangential component, and polynomials
                                         // of degree P_INIT + 1 are present in element interiors. P_INIT = 0
                                         // is for Whitney elements.
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
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
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1.0;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows
                                         // over this limit. This is to prevent h-adaptivity to go on forever.

// Problem parameters.
const double mu_r   = 1.0;
const double kappa  = 1.0;
const double lambda = 1.0;

// Bessel functions, exact solution, and weak forms.
#include "forms.cpp"

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == 1 || marker == 6)
    return BC_ESSENTIAL; // perfect conductor
  else
    return BC_NATURAL;   // impedance
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

int main(int argc, char* argv[])
{
  // Time measurement
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("lshape3q.mesh", &mesh);    // quadrilaterals
  //mloader.load("lshape3t.mesh", &mesh);  // triangles

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Create an Hcurl space with default shapeset.
  HcurlSpace space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form), H2D_SYM);
  wf.add_matrix_form_surf(callback(bilinear_form_surf));
  wf.add_vector_form_surf(linear_form_surf, linear_form_surf_ord);

  // Initialize views.
  OrderView  ordview("Coarse mesh", 600, 0, 600, 500);
  VectorView vecview("Electric Field - VectorView", 0, 0, 600, 500);

  /*
  // View the basis functions.
  VectorBaseView bview;
  vbview.show(&space);
  View::wait(H2DV_WAIT_KEYPRESS);
  */

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof_est, graph_dof_exact, graph_cpu_est, graph_cpu_exact;

  // Initialize refinement selector.
  HcurlProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize the coarse mesh problem.
  LinSystem ls(&wf, &space);

  // Adaptivity loop:
  int as = 1; bool done = false;
  Solution sln_coarse, sln_fine;
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
      int proj_type = 2;    // Hcurl projection.
      ls.project_global(&sln_fine, &sln_coarse, proj_type);
    }

    // Time measurement.
    cpu_time.tick();

    // Calculate error wrt. exact solution.
    info("Calculating error (exact).");
    ExactSolution ex(&mesh, exact);
    double err_exact = hcurl_error(&sln_coarse, &ex) * 100;

    // Show real part of the solution and mesh.
    ordview.show(&space);
    RealFilter real(&sln_coarse);
    vecview.set_min_max_range(0, 1);
    vecview.show(&real, H2D_EPS_HIGH);

    // Skip exact error calculation and visualization time. 
    cpu_time.tick(H2D_SKIP);

    // Calculate error estimate wrt. fine mesh solution.
    info("Calculating error (est).");
    HcurlAdapt hp(&ls);
    hp.set_solutions(&sln_coarse, &sln_fine);
    double err_est_adapt = hp.calc_error() * 100;
    double err_est_hcurl = hcurl_error(&sln_coarse, &sln_fine) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est: %g%%, err_exact: %g%%", 
         ls.get_num_dofs(), rs.get_num_dofs(), err_est_hcurl, err_exact);

    // Add entries to DOF convergence graphs.
    graph_dof_exact.add_values(ls.get_num_dofs(), err_exact);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_dof_est.add_values(ls.get_num_dofs(), err_est_hcurl);
    graph_dof_est.save("conv_dof_est.dat");

    // Add entries to CPU convergence graphs.
    graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact);
    graph_cpu_exact.save("conv_cpu_exact.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_hcurl);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est_adapt too large, adapt the mesh.
    if (err_est_adapt < ERR_STOP) done = true;
    else {
      info("Adapting coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      if (ls.get_num_dofs() >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (!done);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the fine mesh solution - the final result
  vecview.set_title("Final solution");
  vecview.show(&sln_fine);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

