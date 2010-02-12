#include "hermes2d.h"
#include "solver_umfpack.h"

// This example is a continuation of example 08-time-dep. This time, the
// time-dependent laminar incompressible Navier-Stokes equations are
// discretized in time via the implicit Euler method and the Newton's
// iteration is performed in every time step. We encourage you to compare
// the naive linearization method with the Newton's iteration via the
// NEWTON_ON option.
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
// The following parameters can be changed:
//

int NEWTON_ON = 1;             // if NEWTON_ON == 1 then the Newton's iteration is performed
                               // in every time step. Otherwise the convective term is linearized
                               // using the velocities from the previous time step
double RE = 200.0;             // Reynolds number
double VEL_INLET = 1.0;        // inlet velocity (reached after STARTUP_TIME)
double STARTUP_TIME = 1.0;     // during this time, inlet velocity increases gradually
                               // from 0 to VEL_INLET, then it stays constant
double TAU = 0.1;              // time step
double FINAL_TIME = 30000.0;   // length of time interval
int P_INIT_VEL = 2;            // initial polynomial degree for velocity components
int P_INIT_PRESSURE = 1;       // initial polynomial degree for pressure
                               // Note: P_INIT_VEL should always be greater than
                               // P_INIT_PRESSURE because of the inf-sup condition
double NEWTON_TOL = 1e-3;      // convergence criterion for the Newton's method
double H = 5;                  // domain height (necessary to define the parabolic
                               // velocity profile at inlet)

// to better understand boundary conditions
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

scalar xvel_bc_value(int marker, double x, double y) {
  if (marker == marker_left) {
    // time-dependent inlet velocity
    //double val_y = VEL_INLET; //constant profile
    double val_y = VEL_INLET * y*(H-y) / (H/2.)/(H/2.); //parabolic profile with peak VEL_INLET at y = H/2
    if (TIME <= STARTUP_TIME) return val_y * TIME/STARTUP_TIME;
    else return val_y;
  }
  else return 0;
}

int yvel_bc_type(int marker) {
  if (marker == marker_right) return BC_NONE;
  else return BC_ESSENTIAL;
}

int press_bc_type(int marker)
  { return BC_NONE; }

////////////////////////////////////////////////////////////////////////////////////////////////////
// bilinear and linear forms corresponding to simple linearization
// of convective term
template<typename Real, typename Scalar>
Scalar bilinear_form_sym_0_0_1_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) / RE + int_u_v<Real, Scalar>(n, wt, u, v) / TAU;
}

template<typename Real, typename Scalar>
Scalar simple_bilinear_form_unsym_0_0_1_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_w_nabla_u_v<Real, Scalar>(n, wt, ext->fn[0], ext->fn[1], u, v);
}

template<typename Real, typename Scalar>
Scalar simple_linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// bilinear and linear forms corresponding to the Newton's method
template<typename Real, typename Scalar>
Scalar newton_bilinear_form_unsym_0_0(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* vx = ext->fn[0];
  Func<Scalar>* vy = ext->fn[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((vx->val[i] * u->dx[i] + vy->val[i] * u->dy[i]) * v->val[i] + u->val[i] * v->val[i] * vx->dx[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_unsym_0_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* vx = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] * vx->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_unsym_1_0(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* vy = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] * vy->dx[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_unsym_1_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* vx = ext->fn[0];
  Func<Scalar>* vy = ext->fn[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((vx->val[i] * u->dx[i] + vy->val[i] * u->dy[i]) * v->val[i] + u->val[i] * v->val[i] * vy->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_F_0(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* xpr = ext->fn[0];  Func<Scalar>* ypr = ext->fn[1];
  Func<Scalar>* xit = ext->fn[2];  Func<Scalar>* yit = ext->fn[3];  Func<Scalar>* pit = ext->fn[4];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((xit->val[i] - xpr->val[i]) * v->val[i] / TAU +
                       (xit->dx[i] * v->dx[i] + xit->dy[i] * v->dy[i]) / RE +
                       (xit->val[i] * xit->dx[i] + yit->val[i] * xit->dy[i]) * v->val[i] -
                       (pit->val[i] * v->dx[i]));
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_F_1(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* xpr = ext->fn[0];  Func<Scalar>* ypr = ext->fn[1];
  Func<Scalar>* xit = ext->fn[2];  Func<Scalar>* yit = ext->fn[3];  Func<Scalar>* pit = ext->fn[4];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((yit->val[i] - ypr->val[i]) * v->val[i] / TAU +
                       (yit->dx[i] * v->dx[i] + yit->dy[i] * v->dy[i]) / RE +
                       (xit->val[i] * yit->dx[i] + yit->val[i] * yit->dy[i]) * v->val[i] -
                       (pit->val[i] * v->dy[i]));
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_F_2(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* xit = ext->fn[0];  Func<Scalar>* yit = ext->fn[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * (xit->dx[i] * v->val[i] + yit->dy[i] * v->val[i]);
  return result;
}


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // a-priori mesh refinements
  mesh.refine_all_elements();
  mesh.refine_towards_boundary(5, 4, false);
  mesh.refine_towards_boundary(1, 4);
  mesh.refine_towards_boundary(3, 4);

  // initialize the shapeset and the cache
  H1ShapesetBeuchler shapeset;
  PrecalcShapeset pss(&shapeset);

  // spaces for velocities and pressure
  H1Space xvel(&mesh, &shapeset);
  H1Space yvel(&mesh, &shapeset);
  H1Space press(&mesh, &shapeset);

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

  // velocities from the previous time step and for the Newton's iteration
  Solution xprev, yprev, xiter, yiter, piter;
  xprev.set_zero(&mesh);
  yprev.set_zero(&mesh);
  xiter.set_zero(&mesh);
  yiter.set_zero(&mesh);
  piter.set_zero(&mesh);

  // set up weak formulation
  WeakForm wf(3);
  if (NEWTON_ON) {
    wf.add_biform(0, 0, callback(bilinear_form_sym_0_0_1_1), SYM);
    wf.add_biform(0, 0, callback(newton_bilinear_form_unsym_0_0), UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(0, 1, callback(newton_bilinear_form_unsym_0_1), UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(0, 2, callback(bilinear_form_unsym_0_2), ANTISYM);
    wf.add_biform(1, 0, callback(newton_bilinear_form_unsym_1_0), UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(1, 1, callback(bilinear_form_sym_0_0_1_1), SYM);
    wf.add_biform(1, 1, callback(newton_bilinear_form_unsym_1_1), UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(1, 2, callback(bilinear_form_unsym_1_2), ANTISYM);
    wf.add_liform(0, callback(newton_F_0), ANY, 5, &xprev, &yprev, &xiter, &yiter, &piter);
    wf.add_liform(1, callback(newton_F_1), ANY, 5, &xprev, &yprev, &xiter, &yiter, &piter);
    wf.add_liform(2, callback(newton_F_2), ANY, 2, &xiter, &yiter);
  }
  else {
    wf.add_biform(0, 0, callback(bilinear_form_sym_0_0_1_1), SYM);
    wf.add_biform(0, 0, callback(simple_bilinear_form_unsym_0_0_1_1), UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(1, 1, callback(bilinear_form_sym_0_0_1_1), SYM);
    wf.add_biform(1, 1, callback(simple_bilinear_form_unsym_0_0_1_1), UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(0, 2, callback(bilinear_form_unsym_0_2), ANTISYM);
    wf.add_biform(1, 2, callback(bilinear_form_unsym_1_2), ANTISYM);
    wf.add_liform(0, callback(simple_linear_form), ANY, 1, &xprev);
    wf.add_liform(1, callback(simple_linear_form), ANY, 1, &yprev);
  }

  // visualization
  VectorView vview("velocity [m/s]", 0, 0, 1500, 470);
  ScalarView pview("pressure [Pa]", 0, 530, 1500, 470);
  vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  //pview.set_min_max_range(-0.9, 1.0);
  pview.fix_scale_width(80);
  pview.show_mesh(true);

  // matrix solver
  UmfpackSolver umfpack;

  // linear system
  LinSystem linsys(&wf, &umfpack);

  // nonlinear system
  NonlinSystem nonsys(&wf, &umfpack);

  if (NEWTON_ON) {
    // set up the nonlinear system
    nonsys.set_spaces(3, &xvel, &yvel, &press);
    nonsys.set_pss(1, &pss);
  }
  else {
    // set up the linear system
    linsys.set_spaces(3, &xvel, &yvel, &press);
    linsys.set_pss(1, &pss);
  }

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

    if (NEWTON_ON) {

      Solution xsln, ysln, psln;
      int it = 1;
      double res_l2_norm;
      do
      {
        info("\n---- Time step %d, Newton iter %d ---------------------------------\n", i, it++);

        // assemble and solve
        nonsys.assemble();
        nonsys.solve(3, &xsln, &ysln, &psln);

        res_l2_norm = nonsys.get_residuum_l2_norm();
        info("Residuum L2 norm: %g\n", res_l2_norm);
        // want to see Newtons iterations
        //sprintf(title, "Time level %d, Newton iteration %d", i, it-1);
        //vview.set_title(title);
        //vview.show(&xsln, &ysln, EPS_LOW);
        //pview.show(&psln);
        //pview.wait_for_keypress();

        xiter = xsln;
        yiter = ysln;
        piter = psln;

      }
      while (res_l2_norm > NEWTON_TOL);

      // visualization at the end of the time step
      sprintf(title, "Velocity, time %g", TIME);
      vview.set_title(title);
      vview.show(&xprev, &yprev, EPS_LOW);
      sprintf(title, "Pressure, time %g", TIME);
      pview.set_title(title);
      pview.show(&piter);

      // copying result of the Newton's iteration into xprev, yprev
      xprev.copy(&xiter);
      yprev.copy(&yiter);
    }
    else {
      // assemble and solve
      Solution xsln, ysln, psln;
      linsys.assemble();
      linsys.solve(3, &xsln, &ysln, &psln);

      // visualization
      sprintf(title, "Velocity, time %g", TIME);
      vview.set_title(title);
      vview.show(&xsln, &ysln, EPS_LOW);
      sprintf(title, "Pressure, time %g", TIME);
      pview.set_title(title);
      pview.show(&psln);
      //pview.wait_for_keypress();

      xprev = xsln;
      yprev = ysln;
    }

    // uncomment one of the following lines to generate a series of video frames
    //vview.save_numbered_screenshot("velocity%03d.bmp", i, true);
    //pview.save_numbered_screenshot("pressure%03d.bmp", i, true);
    // the frames can then be converted to a video file with the command
    // mencoder "mf://velocity*.bmp" -mf fps=20 -o velocity.avi -ovc lavc -lavcopts vcodec=mpeg4
  }

  // wait for keyboard or mouse input
  View::wait("Waiting for all views to be closed.");
  return 0;
}
