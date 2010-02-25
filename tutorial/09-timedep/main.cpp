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
  return T_INIT;
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return HEATCAP * RHO * int_u_v<Real, Scalar>(n, wt, u, v) / TAU +
         LAMBDA * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return HEATCAP * RHO * int_u_v<Real, Scalar>(n, wt, ext->fn[0], v) / TAU;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return LAMBDA * ALPHA * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return LAMBDA * ALPHA * temp_ext(TIME) * int_v<Real, Scalar>(n, wt, v);
}


int main(int argc, char* argv[])
{
  // load and refine mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("cathedral.mesh", &mesh);
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

  // set initial condition
  Solution tsln;
  tsln.set_const(&mesh, T_INIT);

  // weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form<double, double>, bilinear_form<Ord, Ord>);
  wf.add_biform_surf(0, 0, bilinear_form_surf<double, double>, bilinear_form_surf<Ord, Ord>, marker_air);
  wf.add_liform(0, linear_form<double, double>, linear_form<Ord, Ord>, ANY, 1, &tsln);
  wf.add_liform_surf(0, linear_form_surf<double, double>, linear_form_surf<Ord, Ord>, marker_air);

  // matrix solver
  UmfpackSolver umfpack;

  // linear system
  LinSystem ls(&wf, &umfpack);
  ls.set_spaces(1, &space);
  ls.set_pss(1, &pss);

  // visualisation
  ScalarView Tview("Temperature", 0, 0, 450, 600);
  char title[100];
  sprintf(title, "Time %3.5f, exterior temperature %3.5f", TIME, temp_ext(TIME));
  Tview.set_min_max_range(0,20);
  Tview.set_title(title);
  Tview.fix_scale_width(3);

  // time stepping
  int nsteps = (int)(FINAL_TIME/TAU + 0.5);
  bool rhsonly = false;
  for(int n = 1; n <= nsteps; n++)
  {

    info("\n---- Time %3.5f, time step %d, ext_temp %g ----------", TIME, n, temp_ext(TIME));

    // assemble and solve
    ls.assemble(rhsonly);
    rhsonly = true;
    ls.solve(1, &tsln);

    // shifting the time variable
    TIME += TAU;

    // visualization of solution
    sprintf(title, "Time %3.2f, exterior temperature %3.5f", TIME, temp_ext(TIME));
    Tview.set_title(title);
    Tview.show(&tsln);
    //Tview.wait_for_keypress();

  }

  // Note: a separate example shows how to create videos

  // wait for a view to be closed
  View::wait();
  return 0;
}
