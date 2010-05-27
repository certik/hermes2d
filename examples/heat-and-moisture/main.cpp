#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

using namespace RefinementSelectors;

// This example solves adaptively a time-dependent coupled problem of heat and moisture 
// transfer in massive concrete walls of a nuclear reactor vessel. 
//
// PDE: Lengthy. See the paper P. Solin, L. Dubcova, J. Kruis: Adaptive hp-FEM with Dynamical 
// Meshes for Transient Heat and Moisture Transfer Problems, J. Comput. Appl. Math. 233 (2010) 3103-3112.
//
// The following parameters can be changed:

const int P_INIT = 1;            // Initial polynomial degrees.
const bool MULTI = true;         // MULTI = true  ... use multi-mesh,
                                 // MULTI = false ... use single-mesh.
                                 // Note: In the single mesh option, the meshes are
                                 // forced to be geometrically the same but the
                                 // polynomial degrees can still vary.
const bool SAME_ORDERS = true;   // SAME_ORDERS = true ... when single-mesh is used,
                                 // this forces the meshes for all components to be
                                 // identical, including the polynomial degrees of
                                 // corresponding elements. When multi-mesh is used,
                                 // this parameter is ignored.
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
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;  // Maximum allowed level of hanging nodes:
                                 // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                 // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                 // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                 // Note that regular meshes are not supported, this is due to
                                 // their notoriously bad performance.
const double CONV_EXP = 1.0;     // Default value is 1.0. This parameter influences the selection of
                                 // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double SPACE_TOL = 0.2;    // Stopping criterion for adaptivity (rel. error tolerance between the
                                 // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;     // Adaptivity process stops when the number of degrees of freedom grows over
                                 // this limit. This is mainly to prevent h-adaptivity to go on forever.

// Time step and simulation time.
const double TAU = 5.*24*60*60;           // time step: 120 hours
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

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
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
  // Load the master mesh.
  Mesh basemesh, temp_mesh, moist_mesh;
  H2DReader mloader;
  mloader.load("domain2.mesh", &basemesh);

  // Create temperature and moisture meshes.
  // This also initializes the multimesh hp-FEM.
  temp_mesh.copy(&basemesh);
  moist_mesh.copy(&basemesh);

  // Initialize the shapeset and the cache.
  H1Shapeset shapeset;
  PrecalcShapeset pss1(&shapeset);
  PrecalcShapeset pss2(&shapeset);

  // Create the temperature space.
  H1Space temp(&temp_mesh, &shapeset);
  temp.set_bc_types(temp_bc_type);
  temp.set_essential_bc_values(essential_bc_values);
  temp.set_uniform_order(P_INIT);

  // Create the moisture space.
  H1Space moist(MULTI ? &moist_mesh : &temp_mesh, &shapeset);
  moist.set_bc_types(moist_bc_type);
  moist.set_uniform_order(P_INIT);

  // Define initial conditions.
  Solution temp_prev, moist_prev;
  temp_prev.set_const(&temp_mesh, TEMP_INITIAL);
  moist_prev.set_const(&moist_mesh, MOIST_INITIAL);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_biform(0, 0, callback(bilinear_form_sym_0_0));
  wf.add_biform(0, 1, callback(bilinear_form_sym_0_1));
  wf.add_biform(1, 1, callback(bilinear_form_sym_1_1));
  wf.add_biform(1, 0, callback(bilinear_form_sym_1_0));
  wf.add_liform(0, callback(linear_form_0), H2D_ANY, 1, &temp_prev);
  wf.add_liform(1, callback(linear_form_1), H2D_ANY, 1, &moist_prev);
  wf.add_biform_surf(0, 0, callback(bilinear_form_surf_0_0_ext), MARKER_EXTERIOR_WALL);
  wf.add_biform_surf(1, 1, callback(bilinear_form_surf_1_1_ext), MARKER_EXTERIOR_WALL);
  wf.add_liform_surf(0, callback(linear_form_surf_0_ext), MARKER_EXTERIOR_WALL);
  wf.add_liform_surf(1, callback(linear_form_surf_1_ext), MARKER_EXTERIOR_WALL);

  // Visualize solution and meshes.
  ScalarView temp_view("Temperature [K]", 0, 0, 500, 700);
  OrderView temp_ord("Temperature mesh", 520, 0, 500, 700);
  ScalarView moist_view("Moisture [-]", 1040, 0, 500, 700);
  OrderView moist_ord("Moisture mesh", 1560, 0, 500, 700);
  temp_view.set_min_max_range(TEMP_INITIAL, TEMP_REACTOR_MAX);
  moist_view.set_min_max_range(MOIST_EXTERIOR, MOIST_INITIAL);
  temp_view.show_mesh(false);
  moist_view.show_mesh(false);

  // Initialize the coarse mesh problem.
  UmfpackSolver umfpack;
  LinSystem ls(&wf, &umfpack);
  ls.set_spaces(2, &temp, &moist);
  ls.set_pss(2, &pss1, &pss2);

  // Solutions.
  Solution temp_sln, moist_sln;
  Solution temp_rsln, moist_rsln;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER, &shapeset);

  // Time stepping loop:
  double comp_time = 0.0;
  static int ts = 1;
  while (CURRENT_TIME < SIMULATION_TIME)
  {
    info("---- Physical time = %g s (%d h, %d d, %d y)",
        (CURRENT_TIME + TAU), (int) (CURRENT_TIME + TAU) / 3600,
        (int) (CURRENT_TIME + TAU) / (3600*24), (int) (CURRENT_TIME + TAU) / (3600*24*364));

    // Uniform mesh derefinement.
    temp_mesh.copy(&basemesh);
    moist_mesh.copy(&basemesh);
    temp.set_uniform_order(P_INIT);
    moist.set_uniform_order(P_INIT);

    // Adaptivity loop (in space):
    int as = 1; bool done = false;
    do
    {
      info("---- Time step %d, adaptivity step %d:", ts, as);

      // Enumerate degrees of freedom and update time-dependent Dirichlet BCs.
      int ndof = assign_dofs(2, &temp, &moist);

      // Solve the coarse mesh problem.
      ls.assemble();
      ls.solve(2, &temp_sln, &moist_sln);

      // Solve the fine mesh problem.
      RefSystem rs(&ls);
      rs.assemble();
      rs.solve(2, &temp_rsln, &moist_rsln);

      // Visualize the solution and meshes.
      char title[100];
      sprintf(title, "Temperature mesh, time = %g days", CURRENT_TIME/86400.);
      temp_ord.set_title(title);
      temp_ord.show(&temp);
      sprintf(title, "Moisture mesh, time = %g days", CURRENT_TIME/86400.);
      moist_ord.set_title(title);
      moist_ord.show(&moist);
      sprintf(title, "Temperature, time = %g days", CURRENT_TIME/86400.);
      temp_view.set_title(title);
      temp_view.show(&temp_sln, H2D_EPS_HIGH);
      sprintf(title, "Moisture, time = %g days", CURRENT_TIME/86400.);
      moist_view.set_title(title);
      moist_view.show(&moist_sln, H2D_EPS_HIGH);

      // Calculate element errors and total error estimate.
      H1Adapt hp(Tuple<Space*>(&temp, &moist));
      hp.set_solutions(Tuple<Solution*>(&temp_sln, &moist_sln), Tuple<Solution*>(&temp_rsln, &moist_rsln));
      hp.set_biform(0, 0, callback(bilinear_form_sym_0_0));
      hp.set_biform(0, 1, callback(bilinear_form_sym_0_1));
      hp.set_biform(1, 0, callback(bilinear_form_sym_1_0));
      hp.set_biform(1, 1, callback(bilinear_form_sym_1_1));
      double space_err_est = hp.calc_error(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;

      // Report results.
      info("ndof_temp_coarse: %d, ndof_temp_fine: %d", 
           temp.get_num_dofs(), rs.get_space(0)->get_num_dofs());
      info("ndof_moist_coarse: %d, ndof_moist_fine: %d", 
           moist.get_num_dofs(), rs.get_space(1)->get_num_dofs());
      info("space_err_est: %g%%", space_err_est);

      // If err_est too large, adapt the mesh.
      if (space_err_est > SPACE_TOL) {
        done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY, SAME_ORDERS);
        ndof = assign_dofs(2, &temp, &moist);
        if (ndof >= NDOF_STOP) done = true;
      }
      else done = true;

      as++;
    }
    while (!done);

    // Update time.
    CURRENT_TIME += TAU;

    // Save solutions for the next time step.
    temp_prev = temp_rsln;
    moist_prev = moist_rsln;

    ts++;
  }

  // wait for all views to be closed
  View::wait();
  return 0;
}
