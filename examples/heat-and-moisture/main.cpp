#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"


// This example solves adaptively a time-dependent coupled problem of heat and moisture 
// transfer in massive concrete walls of a nuclear reactor vessel. 
//
// PDE: Lengthy. See the paper P. Solin, L. Dubcova, J. Kruis: Adaptive hp-FEM with Dynamical 
// Meshes for Transient Heat and Moisture Transfer Problems, J. Comput. Appl. Math. 233 (2010) 3103-3112
//

const int P_INIT = 1;            // Initial polynomial degrees 
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
const int ADAPT_TYPE = 0;        // Type of automatic adaptivity:
                                 // ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                 // ADAPT_TYPE = 1 ... adaptive h-FEM,
                                 // ADAPT_TYPE = 2 ... adaptive p-FEM.
const bool ISO_ONLY = false;     // Isotropic refinement flag (concerns quadrilateral elements only).
                                 // ISO_ONLY = false ... anisotropic refinement of quad elements
                                 // is allowed (default),
                                 // ISO_ONLY = true ... only isotropic refinements of quad elements
                                 // are allowed.
const int MESH_REGULARITY = -1;  // Maximum allowed level of hanging nodes:
                                 // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                 // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                 // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                 // Note that regular meshes are not supported, this is due to
                                 // their notoriously bad performance.
const double CONV_EXP = 1.0;     // Default value is 1.0. This parameter influences the selection of
                                 // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int MAXIMUM_ORDER = 10;    // Maximum allowed element degree
const double SPACE_TOL = 0.1;    // Stopping criterion for adaptivity (rel. error tolerance between the
                                 // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;     // Adaptivity process stops when the number of degrees of freedom grows over
                                 // this limit. This is mainly to prevent h-adaptivity to go on forever.

// time step and simulation time
const double SIMULATION_TIME = 3600*24*365*30;  // (seconds) physical time

// equation parameters
const double c_TT = 2.18e+6;
const double d_TT = 2.1;
const double d_Tw = 2.37e-2;
const double k_TT = 25;
const double c_ww = 24.9;
const double d_wT = 1.78e-10;
const double d_ww = 3.02e-8;
const double k_ww = 1.84e-7;

// initial and boundary conditions
const double TEMP_INITIAL = 293.0;       // (Kelvins)
const double MOIST_INITIAL = 0.9;        // (dimensionless)
const double TEMP_EXTERIOR = 293.0;      // (Kelvins)
const double MOIST_EXTERIOR = 0.55;      // (dimensionless)
const double TEMP_REACTOR_MAX = 550.0;   // (Kelvins)
const double REACTOR_START_TIME = 3600*24;  // (seconds) how long does the reactor
                                         // need to warm up linearly from TEMP_INITIAL
                                         // to TEMP_REACTOR_MAX
// material and boundary markers
const int MARKER_SYMMETRY = 1;           // NOTE: this must be compatible with the mesh file!
const int MARKER_REACTOR_WALL = 2;       // NOTE: this must be compatible with the mesh file!
const int MARKER_EXTERIOR_WALL = 5;      // NOTE: this must be compatible with the mesh file!

// for internal use
double current_time = 0.0;      // (seconds) current physical time
double tau = 24*60*60;           // time step: 24 hours


//// boundary conditions ///////////////////////////////////////////////////////////////////////////

int temp_bc_type(int marker)
  { return (marker == MARKER_REACTOR_WALL) ? BC_ESSENTIAL : BC_NATURAL; }

int moist_bc_type(int marker)
  { return BC_NATURAL; }

scalar temp_bc_value(int marker, double x, double y)
{
  if (marker == MARKER_REACTOR_WALL)
  {
    double current_reactor_temperature = TEMP_REACTOR_MAX;
    if (current_time < REACTOR_START_TIME) {
      current_reactor_temperature = TEMP_INITIAL +
        (current_time/REACTOR_START_TIME)*(TEMP_REACTOR_MAX - TEMP_INITIAL);
    }
    return current_reactor_temperature;
  }
  else return 0;
}


//// weak formulation //////////////////////////////////////////////////////////////////////////////
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh basemesh, temp_mesh, moist_mesh;
  H2DReader mloader;
  mloader.load("domain2.mesh", &basemesh); // master mesh
  temp_mesh.copy(&basemesh);
  moist_mesh.copy(&basemesh);

  // shapeset & pss
  H1ShapesetBeuchler shapeset;
  PrecalcShapeset pss1(&shapeset);
  PrecalcShapeset pss2(&shapeset);

  // space for temperature
  H1Space temp(&temp_mesh, &shapeset);
  temp.set_bc_types(temp_bc_type);
  temp.set_bc_values(temp_bc_value);
  temp.set_uniform_order(P_INIT);

  // space for moisture
  H1Space moist(MULTI ? &moist_mesh : &temp_mesh, &shapeset);
  moist.set_bc_types(moist_bc_type);
  moist.set_uniform_order(P_INIT);

  // initial conditions
  Solution temp_prev, moist_prev;
  temp_prev.set_const(&temp_mesh, TEMP_INITIAL);
  moist_prev.set_const(&moist_mesh, MOIST_INITIAL);

  // set up weak formulation
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

  // visualization
  ScalarView temp_view("Temperature [K]", 0, 0, 500, 700);
  OrderView temp_ord("Temperature mesh", 520, 0, 500, 700);
  ScalarView moist_view("Moisture [-]", 1040, 0, 500, 700);
  OrderView moist_ord("Moisture mesh", 1560, 0, 500, 700);
  temp_view.set_min_max_range(TEMP_INITIAL, TEMP_REACTOR_MAX);
  moist_view.set_min_max_range(MOIST_EXTERIOR, MOIST_INITIAL);
  temp_view.show_mesh(false);
  moist_view.show_mesh(false);

  // set up the linear system
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(2, &temp, &moist);
  sys.set_pss(2, &pss1, &pss2);


  // solutions
  Solution temp_sln, moist_sln;
  Solution temp_rsln, moist_rsln;

  // main loop
  double comp_time = 0.0;
  static int step = 0;
  while (current_time < SIMULATION_TIME)
  {
    info("*** Time step %d, current time = %g s (%d h, %d d, %d y)",
        step, (current_time + tau), (int) (current_time + tau) / 3600,
        (int) (current_time + tau) / (3600*24), (int) (current_time + tau) / (3600*24*364));

    //  initiate meshes
    temp_mesh.copy(&basemesh);
    moist_mesh.copy(&basemesh);
    temp.set_uniform_order(P_INIT);
    moist.set_uniform_order(P_INIT);

    double space_err;
    bool done = false;
    int it = 0;
    do // adaptivity in space
    {
      info("*** Adaptivity step %d:", it++);

      // enumerate degrees of freedom and update time-dependent Dirichlet BCs
      int ndof = assign_dofs(2, &temp, &moist);
      if (ndof >= NDOF_STOP) {
        done = true;
        break;
      }

      // assemble and solve
      sys.assemble();
      sys.solve(2, &temp_sln, &moist_sln);

      // visualisation
      char title[100];
      sprintf(title, "Temperature mesh, time = %g days", current_time/86400.);
      temp_ord.set_title(title);
      temp_ord.show(&temp);
      sprintf(title, "Moisture mesh, time = %g days", current_time/86400.);
      moist_ord.set_title(title);
      moist_ord.show(&moist);
      sprintf(title, "Temperature, time = %g days", current_time/86400.);
      temp_view.set_title(title);
      temp_view.show(&temp_sln, H2D_EPS_HIGH);
      sprintf(title, "Moisture, time = %g days", current_time/86400.);
      moist_view.set_title(title);
      moist_view.show(&moist_sln, H2D_EPS_HIGH);

      // solve the fine (reference) problem
      RefSystem rs(&sys);
      rs.assemble();
      rs.solve(2, &temp_rsln, &moist_rsln);

      // calculate errors and adapt the solution
      H1OrthoHP hp(2, &temp, &moist);
      hp.set_biform(0, 0, callback(bilinear_form_sym_0_0));
      hp.set_biform(0, 1, callback(bilinear_form_sym_0_1));
      hp.set_biform(1, 0, callback(bilinear_form_sym_1_0));
      hp.set_biform(1, 1, callback(bilinear_form_sym_1_1));
      double space_err = hp.calc_error_2(&temp_sln, &moist_sln, &temp_rsln, &moist_rsln) * 100;
      info("Energy error est %g%%", space_err);
      if (space_err > SPACE_TOL) hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE, ISO_ONLY, MESH_REGULARITY, CONV_EXP, MAXIMUM_ORDER, SAME_ORDERS);
      else done = true;
    }
    while (!done);

    // update previous solutions
    current_time += tau; step++;
    temp_prev = temp_rsln;
    moist_prev = moist_rsln;
  }

  View::wait();
  return 0;
}
