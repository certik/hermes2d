#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This example solves adaptively a time-dependent coupled problem of heat and moisture 
// transfer in massive concrete walls of a nuclear reactor vessel (simplified axi-symmetric 
// geometry). 
//
// PDE: Lengthy. See the paper P. Solin, L. Dubcova, J. Kruis: Adaptive hp-FEM with Dynamical 
// Meshes for Transient Heat and Moisture Transfer Problems, J. Comput. Appl. Math. 233 (2010) 3103-3112.
//
// The following parameters can be changed:

const bool SOLVE_ON_COARSE_MESH = false; // If true, coarse mesh FE problem is solved in every adaptivity step.
                                         // If false, projection of the fine mesh solution on the coarse mesh is used. 
const int P_INIT = 1;                    // Initial polynomial degrees.
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
const double SPACE_ERR_STOP = 0.5;       // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit. This is mainly to prevent h-adaptivity to go on forever.

// Time step and simulation time.
const double TAU = 5.*24*60*60;                 // time step: 120 hours
const double SIMULATION_TIME = 3600*24*365*30;  // (seconds) physical time

// Equation parameters.
const double c_TT = 2.18e+6;
const double d_TT = 2.1;
const double d_Tw = 2.37e-2;
const double k_TT = 25;
const double c_ww = 24.9;
const double d_wT = 1.78e-10;
const double d_ww = 3.02e-8;
const double k_ww = 1.84e-7;

// Initial and boundary conditions.
const double TEMP_INITIAL = 293.0;           // (Kelvins)
const double MOIST_INITIAL = 0.9;            // (dimensionless)
const double TEMP_EXTERIOR = 293.0;          // (Kelvins)
const double MOIST_EXTERIOR = 0.55;          // (dimensionless)
const double TEMP_REACTOR_MAX = 550.0;       // (Kelvins)
const double REACTOR_START_TIME = 3600*24;   // (seconds) how long does the reactor
                                             // need to warm up linearly from TEMP_INITIAL
                                             // to TEMP_REACTOR_MAX
// Materials and boundary markers.
const int MARKER_SYMMETRY = 1;               
const int MARKER_REACTOR_WALL = 2;           
const int MARKER_EXTERIOR_WALL = 5;          

// Physical time in seconds.
double CURRENT_TIME = 0.0;

// Boundary condition types.
BCType temp_bc_type(int marker)
  { return (marker == MARKER_REACTOR_WALL) ? BC_ESSENTIAL : BC_NATURAL; }

BCType moist_bc_type(int marker)
  { return BC_NATURAL; }

// Essential (Dirichlet) boundary condition values for T.
scalar essential_bc_values_T(int ess_bdy_marker, double x, double y)
{
  if (ess_bdy_marker == MARKER_REACTOR_WALL)
  {
    double current_reactor_temperature = TEMP_REACTOR_MAX;
    if (CURRENT_TIME < REACTOR_START_TIME) {
      current_reactor_temperature = TEMP_INITIAL +
        (CURRENT_TIME/REACTOR_START_TIME)*(TEMP_REACTOR_MAX - TEMP_INITIAL);
    }
    return current_reactor_temperature;
  }
  else return 0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh basemesh, mesh_T, mesh_M;
  H2DReader mloader;
  mloader.load("domain2.mesh", &basemesh);

  // Create temperature and moisture meshes.
  // This also initializes the multimesh hp-FEM.
  mesh_T.copy(&basemesh);
  mesh_M.copy(&basemesh);

  // Create H1 spaces with default shapesets.
  H1Space space_T(&mesh_T, temp_bc_type, essential_bc_values_T, P_INIT);
  H1Space space_M(MULTI ? &mesh_M : &mesh_T, moist_bc_type, NULL, P_INIT);

  // Define constant initial conditions.
  info("Setting initial conditions.");
  Solution T_prev, M_prev;
  T_prev.set_const(&mesh_T, TEMP_INITIAL);
  M_prev.set_const(&mesh_M, MOIST_INITIAL);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(bilinear_form_sym_0_0));
  wf.add_matrix_form(0, 1, callback(bilinear_form_sym_0_1));
  wf.add_matrix_form(1, 1, callback(bilinear_form_sym_1_1));
  wf.add_matrix_form(1, 0, callback(bilinear_form_sym_1_0));
  wf.add_vector_form(0, callback(linear_form_0), H2D_ANY, &T_prev);
  wf.add_vector_form(1, callback(linear_form_1), H2D_ANY, &M_prev);
  wf.add_matrix_form_surf(0, 0, callback(bilinear_form_surf_0_0_ext), MARKER_EXTERIOR_WALL);
  wf.add_matrix_form_surf(1, 1, callback(bilinear_form_surf_1_1_ext), MARKER_EXTERIOR_WALL);
  wf.add_vector_form_surf(0, callback(linear_form_surf_0_ext), MARKER_EXTERIOR_WALL);
  wf.add_vector_form_surf(1, callback(linear_form_surf_1_ext), MARKER_EXTERIOR_WALL);

  // Initialize views.
  ScalarView temp_view("Temperature [K]", 0, 0, 280, 400);
  OrderView temp_ord("Temperature mesh", 300, 0, 280, 400);
  ScalarView moist_view("Moisture [-]", 600, 0, 280, 400);
  OrderView moist_ord("Moisture mesh", 900, 0, 280, 400);
  temp_view.set_min_max_range(TEMP_INITIAL, TEMP_REACTOR_MAX);
  moist_view.set_min_max_range(MOIST_EXTERIOR, MOIST_INITIAL);
  temp_view.show_mesh(false);
  moist_view.show_mesh(false);

  // Error estimate and discrete problem size as a function of physical time.
  SimpleGraph graph_time_err, graph_time_dof;

  // Initialize the coarse mesh problem.
  LinSystem ls(&wf, Tuple<Space*>(&space_T, &space_M));

  // Solutions.
  Solution T_coarse, M_coarse, T_fine, M_fine;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Time stepping loop:
  double comp_time = 0.0;
  static int ts = 1;
  while (CURRENT_TIME < SIMULATION_TIME)
  {
    info("Physical time = %g s (%d h, %d d, %d y)",
        (CURRENT_TIME + TAU), (int) (CURRENT_TIME + TAU) / 3600,
        (int) (CURRENT_TIME + TAU) / (3600*24), (int) (CURRENT_TIME + TAU) / (3600*24*364));

    // Uniform mesh derefinement.
    mesh_T.copy(&basemesh);
    mesh_M.copy(&basemesh);
    space_T.set_uniform_order(P_INIT);
    space_M.set_uniform_order(P_INIT);

    // Adaptivity loop (in space):
    bool done = false;
    double space_err_est;
    int as = 1;
    do
    {
      info("---- Time step %d, adaptivity step %d:", ts, as);

      // Update time-dependent Dirichlet BCs.
      ls.update_essential_bc_values();

      // Solve the fine mesh problem.
      RefSystem rs(&ls);
      info("Solving on fine meshes.");
      rs.assemble();
      rs.solve(Tuple<Solution*>(&T_fine, &M_fine));

      // Either solve on coarse mesh or project the fine mesh solution 
      // on the coarse mesh.
      if (SOLVE_ON_COARSE_MESH) {
        info("Solving on coarse meshes.");
        ls.assemble();
        ls.solve(Tuple<Solution*>(&T_coarse, &M_coarse));
      }
      else {
        info("Projecting fine mesh solutions on coarse meshes.");
        ls.project_global(Tuple<MeshFunction*>(&T_fine, &M_fine), 
                          Tuple<Solution*>(&T_coarse, &M_coarse));
      }

      // Calculate error estimates.
      info("Calculating errors.");
      double T_err_est = h1_error(&T_coarse, &T_fine) * 100;
      double M_err_est = h1_error(&M_coarse, &M_fine) * 100;
      info("T: ndof_coarse: %d, ndof_fine: %d, err_est: %g %%", 
	   ls.get_num_dofs(0), rs.get_num_dofs(0), T_err_est);
      info("M: ndof_coarse: %d, ndof_fine: %d, err_est: %g %%", 
	   ls.get_num_dofs(1), rs.get_num_dofs(1), M_err_est);

      // Calculate errors for adaptivity.
      H1Adapt hp(&ls);
      hp.set_solutions(Tuple<Solution*>(&T_coarse, &M_coarse), 
                       Tuple<Solution*>(&T_fine, &M_fine));
      hp.set_biform(0, 0, callback(bilinear_form_sym_0_0));
      hp.set_biform(0, 1, callback(bilinear_form_sym_0_1));
      hp.set_biform(1, 0, callback(bilinear_form_sym_1_0));
      hp.set_biform(1, 1, callback(bilinear_form_sym_1_1));
      space_err_est = hp.calc_error(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;

      // If err_est too large, adapt the mesh.
      if (space_err_est > SPACE_ERR_STOP) {
        info("Adapting coarse meshes.");
        done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (ls.get_num_dofs() >= NDOF_STOP) done = true;
      }
      else done = true;

      as++;
    }
    while (!done);

    // Visualize the solution and meshes.
    char title[100];
    sprintf(title, "T mesh, time = %g days", CURRENT_TIME/86400.);
    temp_ord.set_title(title);
    temp_ord.show(&space_T);
    sprintf(title, "M mesh, time = %g days", CURRENT_TIME/86400.);
    moist_ord.set_title(title);
    moist_ord.show(&space_M);
    sprintf(title, "T, time = %g days", CURRENT_TIME/86400.);
    temp_view.set_title(title);
    temp_view.show(&T_coarse, H2D_EPS_HIGH);
    sprintf(title, "M, time = %g days", CURRENT_TIME/86400.);
    moist_view.set_title(title);
    moist_view.show(&M_coarse, H2D_EPS_HIGH);

    // Add entries to convergence graphs.
    graph_time_err.add_values(ts*TAU, space_err_est);
    graph_time_err.save("time_error.dat");
    graph_time_dof.add_values(ts*TAU, ls.get_num_dofs());
    graph_time_dof.save("time_dof.dat");

    // Update time.
    CURRENT_TIME += TAU;

    // Save solutions for the next time step.
    T_prev = T_fine;
    M_prev = M_fine;

    ts++;
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
