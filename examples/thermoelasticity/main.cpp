#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

// A massive hollow conductor is heated by induction and cooled by water running inside.
// We will model this problem using linear thermoelasticity equations, where the x-displacement,
// y-displacement, and the temperature will be approximated on individual meshes equipped
// with mutually independent adaptivity mechanisms. Use MULTI = true to use multimesh,
// MULTI = false for single-mesh (all solution components on the samemesh).
//
// PDE: Linear thermoelasticity
//
// BC: u_1 = u_2 = 0 on Gamma_1
//     du_1/dn = du_2/dn = 0 elsewhere
//     temp = TEMP_INNER on Gamma_4
//     negative heat flux with HEAT_FLUX_OUTER elsewhere

const int P_INIT_TEMP = 1;       // Initial polynomial degrees in temperature mesh.
const int P_INIT_DISP = 1;       // Initial polynomial degrees for displacement meshes.
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
const int MAXIMUM_ORDER = 10;     // Maximum allowed element degree
const double ERR_STOP = 0.01;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                 // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;     // Adaptivity process stops when the number of degrees of freedom grows over
                                 // this limit. This is mainly to prevent h-adaptivity to go on forever.

// Problem constants
double HEAT_SRC = 10000.0;       // heat source in the material (caused by induction heating)
double TEMP_INNER = 50;
double HEAT_FLUX_OUTER = -50;
const double E = 2e11;           // steel: E=200 GPa
const double nu = 0.3;
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));
const double l2m = lambda + 2*mu;
const double rho = 8000;
const double g = 9.81;
const double alpha = 13e-6;      // see http://hyperphysics.phy-astr.gsu.edu/hbase/tables/thexp.html

//  Boundary markers:
const int marker_bottom = 1;
const int marker_sides = 2;
const int marker_top = 3;
const int marker_holes = 4;

// Boundary condition types
int bc_types_x(int marker)
  { return (marker == marker_bottom) ? BC_ESSENTIAL : BC_NATURAL; }

int bc_types_y(int marker)
  { return (marker == marker_bottom) ? BC_ESSENTIAL : BC_NATURAL; }

int bc_types_t(int marker)
  { return (marker == marker_holes) ? BC_ESSENTIAL : BC_NATURAL; }

// Boundary condition values
double bc_values_t(int marker, double x, double y)
  { return TEMP_INNER; }

// Weak forms
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh
  Mesh xmesh, ymesh, tmesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &xmesh); // master mesh
  ymesh.copy(&xmesh);                // ydisp will share master mesh with xdisp
  tmesh.copy(&xmesh);                // temp will share master mesh with xdisp

  // Initialize the shapeset and the cache
  H1ShapesetOrtho shapeset;
  PrecalcShapeset xpss(&shapeset);
  PrecalcShapeset ypss(&shapeset);
  PrecalcShapeset tpss(&shapeset);

  // Create the x displacement space
  H1Space xdisp(&xmesh, &shapeset);
  xdisp.set_bc_types(bc_types_x);
  xdisp.set_uniform_order(P_INIT_DISP);

  // Create the y displacement space
  H1Space ydisp(MULTI ? &ymesh : &xmesh, &shapeset);
  ydisp.set_bc_types(bc_types_y);
  ydisp.set_uniform_order(P_INIT_DISP);

  // Create the temperature space
  H1Space temp(MULTI ? &tmesh : &xmesh, &shapeset);
  temp.set_bc_types(bc_types_t);
  temp.set_bc_values(bc_values_t);
  temp.set_uniform_order(P_INIT_TEMP);

  // Initialize the weak formulation
  WeakForm wf(3);
  wf.add_biform(0, 0, callback(bilinear_form_0_0));
  wf.add_biform(0, 1, callback(bilinear_form_0_1), H2D_SYM);
  wf.add_biform(0, 2, callback(bilinear_form_0_2));
  wf.add_biform(1, 1, callback(bilinear_form_1_1));
  wf.add_biform(1, 2, callback(bilinear_form_1_2));
  wf.add_biform(2, 2, callback(bilinear_form_2_2));
  wf.add_liform(1, callback(linear_form_1));
  wf.add_liform(2, callback(linear_form_2));
  wf.add_liform_surf(2, callback(linear_form_surf_2));

  // Visualization
  OrderView xord("X displacement poly degrees", 0, 0, 850, 400);
  OrderView yord("Y displacement poly degrees", 0, 455, 850, 400);
  OrderView tord("Temperature poly degrees", 0, 885, 850, 400);
  ScalarView sview("Von Mises stress [Pa]", 860, 0, 850, 400);
  ScalarView tview("Temperature [deg C]", 860, 455, 850, 400);

  // Matrix solver
  UmfpackSolver solver;

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof, graph_cpu;

  // Adaptivity loop
  int it = 1;
  bool done = false;
  TimePeriod cpu_time;
  Solution x_sln_coarse, y_sln_coarse, t_sln_coarse;
  Solution x_sln_fine, y_sln_fine, t_sln_fine;
  do
  {
    info("!---- Adaptivity step %d ---------------------------------------------", it); it++;

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // enumerate degrees of freedom
    int ndof = assign_dofs(3, &xdisp, &ydisp, &temp);

    // solve the coarse mesh problem
    LinSystem ls(&wf, &solver);
    ls.set_spaces(3, &xdisp, &ydisp, &temp);
    ls.set_pss(3, &xpss, &ypss, &tpss);
    ls.assemble();
    ls.solve(3, &x_sln_coarse, &y_sln_coarse, &t_sln_coarse);

    // time measurement
    cpu_time.tick();

    // report number of dofs
    info("xdof=%d, ydof=%d, tdof=%d", xdisp.get_num_dofs(), ydisp.get_num_dofs(), temp.get_num_dofs());

    // view the solution -- this can be slow; for illustration only
    xord.show(&xdisp);
    yord.show(&ydisp);
    tord.show(&temp);
    VonMisesFilter mises(&x_sln_coarse, &y_sln_coarse, mu, lambda);
    sview.set_min_max_range(0, 4e9);
    sview.show(&mises, H2D_EPS_HIGH);
    tview.show(&t_sln_coarse, H2D_EPS_HIGH);

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // solve the fine mesh problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(3, &x_sln_fine, &y_sln_fine, &t_sln_fine);

    // calculate element errors and total error estimate
    H1OrthoHP hp(3, &xdisp, &ydisp, &temp);
    hp.set_biform(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
    hp.set_biform(0, 1, bilinear_form_0_1<scalar, scalar>, bilinear_form_0_1<Ord, Ord>);
    hp.set_biform(0, 2, bilinear_form_0_2<scalar, scalar>, bilinear_form_0_2<Ord, Ord>);
    hp.set_biform(1, 0, bilinear_form_1_0<scalar, scalar>, bilinear_form_1_0<Ord, Ord>);
    hp.set_biform(1, 1, bilinear_form_1_1<scalar, scalar>, bilinear_form_1_1<Ord, Ord>);
    hp.set_biform(1, 2, bilinear_form_1_2<scalar, scalar>, bilinear_form_1_2<Ord, Ord>);
    hp.set_biform(2, 2, bilinear_form_2_2<scalar, scalar>, bilinear_form_2_2<Ord, Ord>);
    double err_est = hp.calc_error_n(3, &x_sln_coarse, &y_sln_coarse, &t_sln_coarse, &x_sln_fine, &y_sln_fine, &t_sln_fine) * 100;

    // time measurement
    cpu_time.tick();

    // report results
    info("Estimate of error: %g%%", err_est);

    // add entry to DOF convergence graph
    graph_dof.add_values(x_sln_coarse.get_num_dofs() + y_sln_coarse.get_num_dofs() + t_sln_coarse.get_num_dofs(), err_est);
    graph_dof.save("conv_dof.dat");

    // add entry to CPU convergence graph
    graph_cpu.add_values(cpu_time.accumulated(), err_est);
    graph_cpu.save("conv_cpu.dat");

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE, ISO_ONLY, MESH_REGULARITY, CONV_EXP, MAXIMUM_ORDER, SAME_ORDERS);
      ndof = assign_dofs(3, &xdisp, &ydisp, &temp);
      if (ndof >= NDOF_STOP) done = true;
    }

    // time measurement
    cpu_time.tick();
  }
  while (!done);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // show the fine solution - this is the final result
  VonMisesFilter stress_fine(&x_sln_fine, &y_sln_fine, mu, lambda);
  sview.set_title("Final solution");
  sview.set_min_max_range(0, 3e4);
  sview.show(&stress_fine);

  // wait for all views to be closed
  View::wait();
  return 0;
};

