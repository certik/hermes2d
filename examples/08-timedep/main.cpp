#include "hermes2d.h"
#include "solver_umfpack.h"

//  This example shows how to solve a time-dependent PDE discretized
//  in time via the implicit Euler method. The St. Vitus Cathedral
//  in Prague (http://en.wikipedia.org/wiki/St._Vitus_Cathedral)
//  responds to changes in the surrounding air temperature
//  during one 24-hour cycle.
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
const double T_INIT = 10;        // temperature of bottom and initial temperature
const double ALPHA = 10;         // heat flux coefficient for Newton's boundary condition
const double LAMBDA = 50;        // thermal conductivity of the material
const double HEATCAP = 1e6;      // heat capacity
const double RHO = 1;            // material density
const int INIT_REF_NUM = 4;      // number of initial uniform refinements
const double TAU = 300.0;        // time step
const double FINAL_TIME = 86400; // length of time interval (24 hours)
double TIME = 0;

// time-dependent exterior temperature
double temp_ext(double t) {
  return T_INIT + 10. * sin(2*M_PI*t/FINAL_TIME);
}

// boundary markers
int marker_bottom = 1;
int marker_wall = 2;

/********** boundary conditions ***********/

int bc_types(int marker)
{
  if (marker == marker_wall) return BC_NATURAL;
  else return BC_ESSENTIAL;
}

scalar bc_values(int marker, double x, double y)
{
  if (marker == marker_wall) return 0;
  else return T_INIT;
}

/********** linear and bilinear forms ***********/

// previous time step solution
Solution Tprev;

// Implicit Euler method ... dT/dt approximated by (T - T_prev)/TAU
// volumetric forms
scalar bilinear_form_0_0_euler(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return HEATCAP * RHO * int_u_v(fu, fv, ru, rv) / TAU
    + LAMBDA * int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar linear_form_0_euler(RealFunction* fv, RefMap* rv)
{
  return HEATCAP * RHO * int_u_v(&Tprev, fv, Tprev.get_refmap(), rv) / TAU;
}

// surface forms
scalar bilinear_form_0_0_surf(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv, EdgePos *ep)
{
  return LAMBDA * ALPHA * surf_int_u_v(fu, fv, ru, rv, ep);
}

scalar linear_form_0_surf(RealFunction* fv, RefMap* rv, EdgePos *ep)
{
  return LAMBDA * ALPHA * temp_ext(TIME) * surf_int_v(fv, rv, ep);
}

// *************************************************************

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
  wf.add_biform(0, 0, bilinear_form_0_0_euler, UNSYM, ANY, 1, &Tprev);
  wf.add_liform(0, linear_form_0_euler, ANY, 1, &Tprev);
  wf.add_biform_surf(0, 0, bilinear_form_0_0_surf, ANY, 0, marker_wall);
  wf.add_liform_surf(0, linear_form_0_surf, ANY, 0, marker_wall);

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
  Solution T_new;
  int nsteps = (int)(FINAL_TIME/TAU + 0.5);
  for(int n = 1; n <= nsteps; n++)
  {

    info("\n---- Time %3.5f, time step %d, ext_temp %g ----------", TIME, n, temp_ext(TIME));

    // assemble and solve
    ls.assemble();
    ls.solve(1, &T_new);

    // visualization of solution
    sprintf(title, "Time %3.2f, exterior temperature %3.5f", TIME, temp_ext(TIME));
    Tview.set_title(title);
    Tview.show(&T_new);
    //Tview.wait_for_keypress();

    // copying the T_new into T_prev
    Tprev.copy(&T_new);

    TIME += TAU;
  }

  // Note: a separate example shows how to create videos

  printf("Click into the image window and press 'q' to finish.\n");
  View::wait();
  return 0;
}
