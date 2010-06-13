#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

//  This example shows how to solve a time-dependent PDE discretized
//  in time via the implicit Euler method. The St. Vitus Cathedral
//  in Prague (http://en.wikipedia.org/wiki/St._Vitus_Cathedral)
//  responds to changes in the surrounding air temperature
//  during one 24-hour cycle. You will also learn how to use the
//  solution calculated in the previous time step.
//
//  PDE: non-stationary heat transfer equation
//  HEATCAP * RHO * dT/dt - LAMBDA * Laplace T = 0.
//
//  Domain: St. Vitus cathedral (cathedral.mesh).
//
//  IC:  T = T_INIT.
//  BC:  T = T_INIT on the bottom edge ... Dirichlet,
//       dT/dn = ALPHA*(t_exterior(time) - T) ... Newton, time-dependent.
//
//  Time-stepping: implicit Euler.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 4;      // Number of initial uniform mesh refinements.
const int P_INIT = 1;            // Polynomial degree of all mesh elements.
const double TAU = 300.0;        // Time step in seconds.

// Problem parameters.
const double T_INIT = 10;        // Temperature of the ground (also initial temperature).
const double ALPHA = 10;         // Heat flux coefficient for Newton's boundary condition.
const double LAMBDA = 1e5;       // Thermal conductivity of the material.
const double HEATCAP = 1e6;      // Heat capacity.
const double RHO = 3000;         // Material density.
const double FINAL_TIME = 86400; // Length of time interval (24 hours) in seconds.

// Global time variable.
double TIME = 0;

// Time-dependent exterior temperature.
double temp_ext(double t) {
  return T_INIT + 10. * sin(2*M_PI*t/FINAL_TIME);
}

// Boundary markers.
int marker_ground = 1;
int marker_air = 2;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == marker_ground) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return T_INIT;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("cathedral.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(2, 5);

  // Initialize the shapeset.
  H1Shapeset shapeset;

  // Initialize an H1 space.
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_essential_bc_values(essential_bc_values);
  space.set_uniform_order(P_INIT);

  // Set initial condition.
  Solution tsln;
  tsln.set_const(&mesh, T_INIT);

  // Initialize weak formulation.
  WeakForm wf;
  wf.add_biform(bilinear_form<double, double>, bilinear_form<Ord, Ord>);
  wf.add_biform_surf(bilinear_form_surf<double, double>, bilinear_form_surf<Ord, Ord>, marker_air);
  wf.add_liform(linear_form<double, double>, linear_form<Ord, Ord>, H2D_ANY, 1, &tsln);
  wf.add_liform_surf(linear_form_surf<double, double>, linear_form_surf<Ord, Ord>, marker_air);

  // Matrix solver.
  UmfpackSolver solver;

  // Initialize linear system.
  LinSystem ls(&wf, &solver, &space);

  // Initialize views.
  ScalarView Tview("Temperature", 0, 0, 450, 600);
  char title[100];
  sprintf(title, "Time %3.5f, exterior temperature %3.5f", TIME, temp_ext(TIME));
  Tview.set_min_max_range(0,20);
  Tview.set_title(title);
  Tview.fix_scale_width(3);

  // Time stepping:
  int nsteps = (int)(FINAL_TIME/TAU + 0.5);
  bool rhsonly = false;
  for(int ts = 1; ts <= nsteps; ts++)
  {
    info("---- Time step %d, time %3.5f, ext_temp %g", ts, TIME, temp_ext(TIME));

    // Assemble and solve.
    ls.assemble(rhsonly);
    rhsonly = true;
    ls.solve(&tsln);

    // Update the time variable.
    TIME += TAU;

    // Visualize the solution.
    sprintf(title, "Time %3.2f, exterior temperature %3.5f", TIME, temp_ext(TIME));
    Tview.set_title(title);
    Tview.show(&tsln);
  }

  // Wait for the view to be closed.
  View::wait();
  return 0;
}
