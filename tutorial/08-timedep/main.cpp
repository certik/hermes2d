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
//  HEATCAP * RHO * dT/dt - LAMBDA * Laplace T = 0
//
//  Domain: St. Vitus cathedral (cathedral.mesh)
//
//  IC:  T = T_INIT
//  BC:  T = T_INIT on the bottom edge ... Dirichlet
//       dT/dn = ALPHA*(t_exterior(time) - T) ... Newton, time-dependent
//
//  Time-stepping: implicit Euler
//
//  The following parameters can be played with:

const int P_INIT = 1;            // polynomial degree of elements
const int INIT_REF_NUM = 4;      // number of initial uniform refinements
const double TAU = 300.0;        // time step in seconds

// problem constants
const double T_INIT = 10;        // temperature of the ground (also initial temperature)
const double ALPHA = 10;         // heat flux coefficient for Newton's boundary condition
const double LAMBDA = 1e5;       // thermal conductivity of the material
const double HEATCAP = 1e6;      // heat capacity
const double RHO = 3000;         // material density
const double FINAL_TIME = 86400; // length of time interval (24 hours) in seconds

// global variable
double TIME = 0;

// time-dependent exterior temperature
double temp_ext(double t) {
  return T_INIT + 10. * sin(2*M_PI*t/FINAL_TIME);
}

// boundary markers
int marker_ground = 1;
int marker_air = 2;

// boundary condition types
int bc_types(int marker)
{
  if (marker == marker_ground) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// function values for Dirichlet boundary markers
scalar bc_values(int marker, double x, double y)
{
  if (marker == marker_ground) return T_INIT;
}

// previous time step solution
Solution Tprev;

// Implicit Euler method ... dT/dt approximated by (Tnew - Tprev)/TAU
// bilinear forms
scalar bilinear_form_0_0_euler(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return HEATCAP * RHO * int_u_v(fu, fv, ru, rv) / TAU
    + LAMBDA * int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar bilinear_form_0_0_surf(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv, EdgePos *ep)
{
  return LAMBDA * ALPHA * surf_int_u_v(fu, fv, ru, rv, ep);
}

// linear forms
scalar linear_form_0_euler(RealFunction* fv, RefMap* rv)
{
  return HEATCAP * RHO * int_u_v(&Tprev, fv, Tprev.get_refmap(), rv) / TAU;
}

scalar linear_form_0_surf(RealFunction* fv, RefMap* rv, EdgePos *ep)
{
  return LAMBDA * ALPHA * temp_ext(TIME + TAU) * surf_int_v(fv, rv, ep);
}

int main(int argc, char* argv[])
{
  // load and refine mesh
  Mesh mesh;
  mesh.load("cathedral.mesh");
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(2, 5);

  // set up shapeset
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // set up spaces
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // enumerate basis functions
  space.assign_dofs();

  // weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form_0_0_euler, SYM, ANY);
  wf.add_liform(0, linear_form_0_euler, ANY, 1, &Tprev);
  wf.add_biform_surf(0, 0, bilinear_form_0_0_surf, marker_air);
  wf.add_liform_surf(0, linear_form_0_surf, marker_air);

  // matrix solver
  UmfpackSolver umfpack;

  // linear system
  LinSystem ls(&wf, &umfpack);
  ls.set_spaces(1, &space);
  ls.set_pss(1, &pss);

  // set initial condition
  Tprev.set_const(&mesh, T_INIT);

  // visualisation
  ScalarView Tview("Temperature", 0, 0, 450, 600);
  char title[100];
  sprintf(title, "Time %3.5f, exterior temperature %3.5f", TIME, temp_ext(TIME));
  Tview.set_min_max_range(0,20);
  Tview.set_title(title);
  Tview.fix_scale_width(3);

  // time stepping
  Solution Tnew;
  int nsteps = (int)(FINAL_TIME/TAU + 0.5);
  bool rhsonly = false;
  for(int n = 1; n <= nsteps; n++)
  {

    info("\n---- Time %3.5f, time step %d, ext_temp %g ----------", TIME, n, temp_ext(TIME));

    // assemble and solve
    ls.assemble(rhsonly);
    rhsonly = true;
    ls.solve(1, &Tnew);

    // shifting the time variable
    TIME += TAU;

    // visualization of solution
    sprintf(title, "Time %3.2f, exterior temperature %3.5f", TIME, temp_ext(TIME));
    Tview.set_title(title);
    Tview.show(&Tnew);
    //Tview.wait_for_keypress();

    // copying Tnew into Tprev
    Tprev = Tnew;
  }

  // Note: a separate example shows how to create videos

  // wait for keyboard or mouse input
  printf("Waiting for keyboard or mouse input.\n");
  View::wait();
  return 0;
}
