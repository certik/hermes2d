#include "hermes2d.h"
#include "solver_umfpack.h"

// The time-dependent laminar incompressible Navier-Stokes equations are
// discretized in time via the implicit Euler method. The convective term
// is linearized simply by replacing the velocity in front of the nabla
// operator with the velocity from last time step. Velocity is approximated
// using continuous elements, and pressure by means or discontinuous (L2)
// elements. This makes the velocity discretely divergence-free, meaning
// that the integral of div(v) over every element is zero. The problem
// has a steady symmetric solution which is unstable. Note that after some
// time (around t = 100), numerical errors induce oscillations. The
// approximation becomes unsteady and thus diverges from the exact solution.
// Interestingly, this happens even with a completely symmetric mesh.
//
// PDE: incompressible Navier-Stokes equations in the form
// \partial v / \partial t - \Delta v / Re + (v \cdot \nabla) v + \nabla p = 0,
// div v = 0
//
// BC: u_1 is a time-dependent constant and u_2 = 0 on Gamma_4 (inlet)
//     u_1 = u_2 = 0 on Gamma_1 (bottom), Gamma_3 (top) and Gamma_5 (obstacle)
//     "do nothing" on Gamma_2 (outlet)
//
// TODO: Implement Crank-Nicolson so that comparisons with implicit Euler can be made
//
// The following parameters can be played with:

const double RE = 1000.0;            // Reynolds number
const double VEL_INLET = 1.0;        // inlet velocity (reached after STARTUP_TIME)
const double STARTUP_TIME = 1.0;     // during this time, inlet velocity increases gradually
                                     // from 0 to VEL_INLET, then it stays constant
const double TAU = 0.5;              // time step
const double FINAL_TIME = 3000.0;    // length of time interval
const int P_INIT_VEL = 2;            // polynomial degree for velocity components
const int P_INIT_PRESSURE = 1;       // polynomial degree for pressure
                                     // Note: P_INIT_VEL should always be greater than
                                     // P_INIT_PRESSURE because of the inf-sup condition
const double H = 5.0;                // domain height (necessary to define the parabolic
                                     // velocity profile at inlet)

//  boundary markers
int marker_bottom = 1;
int marker_right  = 2;
int marker_top = 3;
int marker_left = 4;
int marker_obstacle = 5;

// global time variable
double TIME = 0;

// definition of boundary conditions
int xvel_bc_type(int marker) {
  if (marker == marker_right) return BC_NONE;
  else return BC_ESSENTIAL;
}

int yvel_bc_type(int marker) {
  if (marker == marker_right) return BC_NONE;
  else return BC_ESSENTIAL;
}

int press_bc_type(int marker)
  { return BC_NONE; }

scalar xvel_bc_value(int marker, double x, double y) {
  if (marker == marker_left) {
    // time-dependent inlet velocity
//     double val_y = VEL_INLET; //constant profile
    double val_y = VEL_INLET * y*(H-y) / (H/2.)/(H/2.); //parabolic profile with peak VEL_INLET at y = H/2
    if (TIME <= STARTUP_TIME) return val_y * TIME/STARTUP_TIME;
    else return val_y;
    return val_y;
  }
  else return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar bilinear_form_sym_0_0_1_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) / RE + int_u_v<Real, Scalar>(n, wt, u, v) / TAU;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_unsym_0_0_1_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_w_nabla_u_v<Real, Scalar>(n, wt, ext->fn[0], ext->fn[1], u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_u_v<Real, Scalar>(n, wt, ext->fn[0], v) / TAU;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_unsym_0_2(int n, double *wt, Func<Real> *p, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return - int_u_dvdx<Real, Scalar>(n, wt, p, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_unsym_1_2(int n, double *wt, Func<Real> *p, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return - int_u_dvdy<Real, Scalar>(n, wt, p, v);
}


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain-quad.mesh", &mesh); // unstructured triangular mesh available in domain-tri.mesh

  // a-priori mesh refinements
  mesh.refine_all_elements();
  mesh.refine_towards_boundary(5, 4, false);
  mesh.refine_towards_boundary(1, 4);
  mesh.refine_towards_boundary(3, 4);

  // display the mesh
  //MeshView mview("Navier-Stokes Example - Mesh", 100, 100, 1100, 400);
  //mview.show(&mesh);
  //mview.wait_for_keypress();

  // initialize the shapesets and the cache
  H1Shapeset shapeset_h1;
  PrecalcShapeset pss_h1(&shapeset_h1);
  L2Shapeset shapeset_l2;
  PrecalcShapeset pss_l2(&shapeset_l2);

  // H1 spaces for velocities and L2 for pressure
  H1Space xvel(&mesh, &shapeset_h1);
  H1Space yvel(&mesh, &shapeset_h1);
  L2Space press(&mesh, &shapeset_l2);

  // initialize boundary conditions
  xvel.set_bc_types(xvel_bc_type);
  xvel.set_bc_values(xvel_bc_value);
  yvel.set_bc_types(yvel_bc_type);
  press.set_bc_types(press_bc_type);

  // set velocity and pressure polynomial degrees
  xvel.set_uniform_order(P_INIT_VEL);
  yvel.set_uniform_order(P_INIT_VEL);
  press.set_uniform_order(P_INIT_PRESSURE);

  // assign degrees of freedom
  int ndofs = 0;
  ndofs += xvel.assign_dofs(ndofs);
  ndofs += yvel.assign_dofs(ndofs);
  ndofs += press.assign_dofs(ndofs);

  // initial condition: xprev and yprev are zero
  Solution xprev, yprev;
  xprev.set_zero(&mesh);
  yprev.set_zero(&mesh);

  // set up weak formulation
  WeakForm wf(3);
  wf.add_biform(0, 0, callback(bilinear_form_sym_0_0_1_1), SYM);
  wf.add_biform(0, 0, callback(bilinear_form_unsym_0_0_1_1), UNSYM, ANY, 2, &xprev, &yprev);
  wf.add_biform(1, 1, callback(bilinear_form_sym_0_0_1_1), SYM);
  wf.add_biform(1, 1, callback(bilinear_form_unsym_0_0_1_1), UNSYM, ANY, 2, &xprev, &yprev);
  wf.add_biform(0, 2, callback(bilinear_form_unsym_0_2), ANTISYM);
  wf.add_biform(1, 2, callback(bilinear_form_unsym_1_2), ANTISYM);
  wf.add_liform(0, callback(linear_form), ANY, 1, &xprev);
  wf.add_liform(1, callback(linear_form), ANY, 1, &yprev);

  // visualization
  VectorView vview("velocity [m/s]", 0, 0, 1500, 470);
  ScalarView pview("pressure [Pa]", 0, 530, 1500, 470);
  vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  pview.show_mesh(false);
  pview.fix_scale_width(80);
  // fixing scale width (for nicer videos). Note: creation of videos is
  // discussed in a separate example
  vview.fix_scale_width(5);
  pview.fix_scale_width(5);

  // set up the linear system
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(3, &xvel, &yvel, &press);
  sys.set_pss(3, &pss_h1, &pss_h1, &pss_l2);

  // main loop
  char title[100];
  int num_time_steps = FINAL_TIME / TAU;
  for (int i = 1; i <= num_time_steps; i++)
  {
    TIME += TAU;

    info("\n---- Time step %d, time = %g -----------------------------------", i, TIME);

    // this is needed to update the time-dependent boundary conditions
    ndofs = 0;
    ndofs += xvel.assign_dofs(ndofs);
    ndofs += yvel.assign_dofs(ndofs);
    ndofs += press.assign_dofs(ndofs);

    // assemble and solve
    Solution xsln, ysln, psln;
    sys.assemble();
    sys.solve(3, &xsln, &ysln, &psln);

    // visualization
    sprintf(title, "Velocity, time %g", TIME);
    vview.set_title(title);
    vview.show(&xprev, &yprev, EPS_LOW);
    sprintf(title, "Pressure, time %g", TIME);
    pview.set_title(title);
    pview.show(&psln);

    xprev = xsln;
    yprev = ysln;
  }

  // wait for keyboard or mouse input
  View::wait("Waiting for all views to be closed.");
  return 0;
}
