#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

//  This example is a very simple flame propagation model (laminar flame,
//  zero flow velocity), and its purpose is to show how the Newton's method
//  is applied to a time-dependent two-equation system.
//
//  PDEs:
//
//  dT/dt - laplace T = omega(T,Y)
//  dY/dt - 1/Le * laplace Y = - omega(T,Y)
//
//  Domain: rectangle with cooled rods
//
//  BC:  T = 1, Y = 0 on the inlet
//       dT/dn = - kappa T on cooled rods
//       dT/dn = 0, dY/dn = 0 elsewhere
//
//  Time-stepping: second order BDF formula

const int P_INIT = 2;                  // Initial polynomial degree
const double TAU = 0.5;                // Time step
const double T_FINAL = 60.0;           // Time interval length
const int PROJ_TYPE = 1;               // For the projection of the initial condition
                                       // on the initial mesh: 1 = H1 projection, 0 = L2 projection
const int INIT_GLOB_REF_NUM = 2;       // Number of initial uniform mesh refinements
const double NEWTON_TOL = 1e-4;        // Stopping criterion for the Newton's method
const int NEWTON_MAX_ITER = 15;        // Maximum allowed number of Newton iterations

// Problem constants
const double Le    = 1.0;
const double alpha = 0.8;
const double beta  = 10.0;
const double kappa = 0.1;
const double x1    = 9.0;

// Boundary conditions
int bc_types(int marker)
  { return (marker == 1) ? BC_ESSENTIAL : BC_NATURAL; }

scalar temp_bc_values(int marker, double x, double y)
  { return (marker == 1) ? 1.0 : 0; }

// Initial conditions
scalar temp_ic(double x, double y, scalar& dx, scalar& dy)
  { return (x <= x1) ? 1.0 : exp(x1 - x); }

scalar conc_ic(double x, double y, scalar& dx, scalar& dy)
  { return (x <= x1) ? 0.0 : 1.0 - exp(Le*(x1 - x)); }

// Weak forms, definition of reaction rate omega
# include "forms.cpp"

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // initial mesh refinements
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create H1 spaces
  H1Space tspace(&mesh, &shapeset);
  H1Space cspace(&mesh, &shapeset);
  tspace.set_bc_types(bc_types);
  tspace.set_bc_values(temp_bc_values);
  cspace.set_bc_types(bc_types);
  tspace.set_uniform_order(P_INIT);
  cspace.set_uniform_order(P_INIT);

  // enumerate degrees of freedom
  int ndof = assign_dofs(2, &tspace, &cspace);

  // solutions for the Newton's iteration and time stepping
  Solution t_prev_time_1, y_prev_time_1, t_prev_time_2, y_prev_time_2, t_prev_newton, y_prev_newton, tsln, csln;

  // setting initial conditions
  t_prev_time_1.set_exact(&mesh, temp_ic); y_prev_time_1.set_exact(&mesh, conc_ic);
  t_prev_time_2.set_exact(&mesh, temp_ic); y_prev_time_2.set_exact(&mesh, conc_ic);
  t_prev_newton.set_exact(&mesh, temp_ic);  y_prev_newton.set_exact(&mesh, conc_ic);

  // defining filters for the reaction rate omega
  DXDYFilter omega(omega_fn, &t_prev_newton, &y_prev_newton);
  DXDYFilter omega_dt(omega_dt_fn, &t_prev_newton, &y_prev_newton);
  DXDYFilter omega_dy(omega_dy_fn, &t_prev_newton, &y_prev_newton);

  // visualization
  ScalarView rview("Reaction rate", 0, 0, 1600, 460);

  // initialize the weak formulation
  WeakForm wf(2);
  wf.add_biform(0, 0, callback(newton_bilinear_form_0_0), H2D_UNSYM, H2D_ANY, 1, &omega_dt);
  wf.add_biform_surf(0, 0, callback(newton_bilinear_form_0_0_surf), 3);
  wf.add_biform(0, 1, callback(newton_bilinear_form_0_1), H2D_UNSYM, H2D_ANY, 1, &omega_dy);
  wf.add_biform(1, 0, callback(newton_bilinear_form_1_0), H2D_UNSYM, H2D_ANY, 1, &omega_dt);
  wf.add_biform(1, 1, callback(newton_bilinear_form_1_1), H2D_UNSYM, H2D_ANY, 1, &omega_dy);
  wf.add_liform(0, callback(newton_linear_form_0), H2D_ANY, 4, &t_prev_newton, &t_prev_time_1, &t_prev_time_2, &omega);
  wf.add_liform_surf(0, callback(newton_linear_form_0_surf), 3, 1, &t_prev_newton);
  wf.add_liform(1, callback(newton_linear_form_1), H2D_ANY, 4, &y_prev_newton, &y_prev_time_1, &y_prev_time_2, &omega);

  // initialize the nonlinear system and solver
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(2, &tspace, &cspace);
  nls.set_pss(1, &pss);
  nls.set_ic(&t_prev_time_1, &y_prev_time_1, &t_prev_newton, &y_prev_newton, PROJ_TYPE);

  // time stepping loop
  double current_time = 0.0;
  int t_step = 0;
  do {
    t_step++;
    info("**** Time step %d, t = %g s:", t_step, current_time);

    // Newton's method
    if (!nls.solve_newton_2(&t_prev_newton, &y_prev_newton, NEWTON_TOL, NEWTON_MAX_ITER,
                       &omega, &omega_dt, &omega_dy)) error("Newton's method did not converge.");

    // visualization
    DXDYFilter omega_view(omega_fn, &t_prev_newton, &y_prev_newton);
    rview.set_min_max_range(0.0,2.0);
    char title[100];
    sprintf(title, "Reaction rate, t = %g", current_time);
    rview.set_title(title);
    rview.show(&omega_view);

    // update current time
    current_time += TAU;

    // store two time levels of previous solutions
    t_prev_time_2.copy(&t_prev_time_1);
    y_prev_time_2.copy(&y_prev_time_1);
    t_prev_time_1.copy(&t_prev_newton);
    y_prev_time_1.copy(&y_prev_newton);
  } while (current_time <= T_FINAL);

  // wait for all views to be closed
  View::wait();
  return 0;
}
