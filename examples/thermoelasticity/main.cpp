#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"


using namespace RefinementSelectors;

// A massive hollow conductor is heated by induction and cooled by water running inside.
// We will model this problem using linear thermoelasticity equations, where the x-displacement,
// y-displacement, and the temperature will be approximated on individual meshes equipped
// with mutually independent adaptivity mechanisms. Use MULTI = true to use multimesh,
// MULTI = false for single-mesh (all solution components on the samemesh).
//
// PDE: Linear thermoelasticity.
//
// BC: u_1 = u_2 = 0 on Gamma_1
//     du_1/dn = du_2/dn = 0 elsewhere
//     temp = TEMP_INNER on Gamma_4
//     negative heat flux with HEAT_FLUX_OUTER elsewhere.

const bool SOLVE_ON_COARSE_MESH = false; // If true, coarse mesh FE problem is solved in every adaptivity step.
                                         // If false, projection of the fine mesh solution on the coarse mesh is used. 
const int P_INIT_TEMP = 1;               // Initial polynomial degrees in temperature mesh.
const int P_INIT_DISP = 1;               // Initial polynomial degrees for displacement meshes.
const bool MULTI = true;                 // MULTI = true  ... use multi-mesh,
                                         // MULTI = false ... use single-mesh.
                                         // Note: In the single mesh option, the meshes are
                                         // forced to be geometrically the same but the
                                         // polynomial degrees can still vary.
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
const double ERR_STOP = 0.01;            // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit. This is mainly to prevent h-adaptivity to go on forever.

// Problem parameters.
double HEAT_SRC = 10000.0;               // Heat source in the material (caused by induction heating).
double TEMP_INNER = 50;
double HEAT_FLUX_OUTER = -50;
const double E = 2e11;                   // Steel: E=200 GPa.
const double nu = 0.3;
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));
const double l2m = lambda + 2*mu;
const double rho = 8000;
const double g = 9.81;
const double alpha = 13e-6;              // See http://hyperphysics.phy-astr.gsu.edu/hbase/tables/thexp.html.

//  Boundary markers.
const int marker_bottom = 1;
const int marker_sides = 2;
const int marker_top = 3;
const int marker_holes = 4;

// Boundary condition types.
BCType bc_types_x(int marker)
  { return (marker == marker_bottom) ? BC_ESSENTIAL : BC_NATURAL; }

BCType bc_types_y(int marker)
  { return (marker == marker_bottom) ? BC_ESSENTIAL : BC_NATURAL; }

BCType bc_types_t(int marker)
  { return (marker == marker_holes) ? BC_ESSENTIAL : BC_NATURAL; }

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values_temp(int ess_bdy_marker, double x, double y)
  { return TEMP_INNER; }

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh xmesh, ymesh, tmesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &xmesh); // Master mesh.

  // Initialize multimesh hp-FEM.
  ymesh.copy(&xmesh);                  // Ydisp will share master mesh with xdisp.
  tmesh.copy(&xmesh);                  // Temp will share master mesh with xdisp.

  // Create H1 spaces with default shapesets.
  H1Space xdisp(&xmesh, bc_types_x, NULL, P_INIT_DISP);
  H1Space ydisp(MULTI ? &ymesh : &xmesh, bc_types_y, NULL, P_INIT_DISP);
  H1Space temp(MULTI ? &tmesh : &xmesh, bc_types_t, essential_bc_values_temp, P_INIT_TEMP);

  // Initialize the weak formulation.
  WeakForm wf(3);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0));
  wf.add_matrix_form(0, 1, callback(bilinear_form_0_1), H2D_SYM);
  wf.add_matrix_form(0, 2, callback(bilinear_form_0_2));
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1));
  wf.add_matrix_form(1, 2, callback(bilinear_form_1_2));
  wf.add_matrix_form(2, 2, callback(bilinear_form_2_2));
  wf.add_vector_form(1, callback(linear_form_1));
  wf.add_vector_form(2, callback(linear_form_2));
  wf.add_vector_form_surf(2, callback(linear_form_surf_2));

  // Initialize views.
  OrderView xord("X displacement mesh", 0, 0, 600, 300);
  OrderView yord("Y displacement mesh", 0, 350, 600, 300);
  OrderView tord("Temperature mesh", 0, 675, 600, 300);
  ScalarView sview("Von Mises stress [Pa]", 610, 0, 600, 300);
  ScalarView tview("Temperature [deg C]", 610, 350, 600, 300);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize the coarse mesh problem.
  LinSystem ls(&wf, Tuple<Space*>(&xdisp, &ydisp, &temp));

  // Adaptivity loop:
  int as = 1; bool done = false;
  Solution x_sln_coarse, y_sln_coarse, t_sln_coarse;
  Solution x_sln_fine, y_sln_fine, t_sln_fine;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Solve the fine mesh problems.
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(Tuple<Solution*>(&x_sln_fine, &y_sln_fine, &t_sln_fine));

    // Either solve on coarse meshes or project the fine mesh solutions
    // on the coarse mesh.
    if (SOLVE_ON_COARSE_MESH) {
      info("Solving on coarse mesh.");
      ls.assemble();
      ls.solve(Tuple<Solution*>(&x_sln_coarse, &y_sln_coarse, &t_sln_coarse));
    }
    else {
      info("Projecting fine mesh solution on coarse mesh.");
      ls.project_global(Tuple<MeshFunction*>(&x_sln_fine, &y_sln_fine, &t_sln_fine), 
                        Tuple<Solution*>(&x_sln_coarse, &y_sln_coarse, &t_sln_coarse));
    }

    // Time measurement.
    cpu_time.tick();

    // View the solutions.
    xord.show(&xdisp);
    yord.show(&ydisp);
    tord.show(&temp);
    VonMisesFilter mises(&x_sln_coarse, &y_sln_coarse, mu, lambda);
    sview.set_min_max_range(0, 4e9);
    sview.show(&mises, H2D_EPS_HIGH);
    tview.show(&t_sln_coarse, H2D_EPS_HIGH);

    // Skip visualization time.
    cpu_time.tick(H2D_SKIP);

    // Calculate element errors and total error estimate.
    H1Adapt hp(&ls);
    hp.set_solutions(Tuple<Solution*>(&x_sln_coarse, &y_sln_coarse, &t_sln_coarse), 
                     Tuple<Solution*>(&x_sln_fine, &y_sln_fine, &t_sln_fine));
    hp.set_biform(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
    hp.set_biform(0, 1, bilinear_form_0_1<scalar, scalar>, bilinear_form_0_1<Ord, Ord>);
    hp.set_biform(0, 2, bilinear_form_0_2<scalar, scalar>, bilinear_form_0_2<Ord, Ord>);
    hp.set_biform(1, 0, bilinear_form_1_0<scalar, scalar>, bilinear_form_1_0<Ord, Ord>);
    hp.set_biform(1, 1, bilinear_form_1_1<scalar, scalar>, bilinear_form_1_1<Ord, Ord>);
    hp.set_biform(1, 2, bilinear_form_1_2<scalar, scalar>, bilinear_form_1_2<Ord, Ord>);
    hp.set_biform(2, 2, bilinear_form_2_2<scalar, scalar>, bilinear_form_2_2<Ord, Ord>);
    double err_est = hp.calc_error(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_ABS) * 100;

    // Report results.
    info("ndof_x_coarse: %d, ndof_y_coarse: %d, ndof_t_coarse: %d", 
      ls.get_num_dofs(0), ls.get_num_dofs(1), ls.get_num_dofs(2));
    info("ndof_x_fine: %d, ndof_y_fine: %d, ndof_t_fine: %d", 
      rs.get_num_dofs(0), rs.get_num_dofs(1), rs.get_num_dofs(2));
    info("err_est: %g%%", err_est);

    // Add entry to DOF convergence graph.
    graph_dof.add_values(ls.get_num_dofs(), err_est);
    graph_dof.save("conv_dof.dat");

    // Add entry to CPU convergence graph.
    graph_cpu.add_values(cpu_time.accumulated(), err_est);
    graph_cpu.save("conv_cpu.dat");

    // If err_est too large, adapt the mesh.
    if (err_est < ERR_STOP) done = true;
    else {
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      if (ls.get_num_dofs() >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (!done);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the fine solution - the final result.
  VonMisesFilter stress_fine(&x_sln_fine, &y_sln_fine, mu, lambda);
  sview.set_title("Final solution");
  sview.show_mesh(false);
  sview.set_min_max_range(0, 3e4);
  sview.show(&stress_fine);

  // Wait for all views to be closed.
  View::wait();
  return 0;
};

