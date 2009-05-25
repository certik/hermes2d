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
                               // using the velocities from the previous time step (as in example
                               // 08-time-dep)
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

// same as int_u_dvdx() but now 'v' is a solution
inline double int_u_dvdx_II(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double *uval = fu->get_fn_values();
  double *dvdx, *dvdy;
  fv->get_dx_dy_values(dvdx, dvdy); // 'v' solution => derivatives already transformed to physical element

  double result = 0.0;
  h1_integrate_dd_expression(uval[i] * dvdx[i]);
  return result;
}

// same as int_u_dvdy() but now 'v' is a solution
inline double int_u_dvdy_II(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double *uval = fu->get_fn_values();
  double *dvdx, *dvdy;
  fv->get_dx_dy_values(dvdx, dvdy); // 'v' solution => derivatives already transformed to physical element

  double result = 0.0;
  h1_integrate_dd_expression(uval[i] * dvdy[i]);
  return result;
}

#define int_dudx_v_II(fu, fv, ru, rv) int_u_dvdx_II(fv, fu, rv, ru)
#define int_dudy_v_II(fu, fv, ru, rv) int_u_dvdy_II(fv, fu, rv, ru)

// same as int_grad_u_grad_v but now fu is a solution
inline double int_grad_u_grad_v_II(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double *dudx, *dudy, *dvdx, *dvdy;
  fu->get_dx_dy_values(dudx, dudy); // 'u' solution => derivatives already transformed to physical element
  fv->get_dx_dy_values(dvdx, dvdy);

  double result = 0.0;
  h1_integrate_dd_expression(dudx[i] * t_dvdx + dudy[i] * t_dvdy);
  return result;
}

// same as int_w_nabla_u_v() but now 'u' is a solution
inline double int_w_nabla_u_v_II(RealFunction* w1, RealFunction* w2,
                              RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = fu->get_fn_order() + fv->get_fn_order() +
          w1->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);

  w1->set_quad_order(o, FN_VAL);
  w2->set_quad_order(o, FN_VAL);
  fu->set_quad_order(o);
  fv->set_quad_order(o, FN_VAL);

  double *dudx, *dudy;
  fu->get_dx_dy_values(dudx, dudy); // 'u' solution => derivatives already transformed to physical element
  double* vval = fv->get_fn_values();
  double* w1val = w1->get_fn_values();
  double* w2val = w2->get_fn_values();

  double result = 0.0;
  h1_integrate_dd_expression((w1val[i] * dudx[i] + w2val[i] * dudy[i]) * vval[i]);
  return result;
}

inline double int_u_dvdx_w(RealFunction* fu, RealFunction* fv, RealFunction* fw, RefMap* ru, RefMap* rv, RefMap* rw)
// assumes that 'v' is a solution
{
  Quad2D* quad = fu->get_quad_2d();

  int o = fu->get_fn_order() + fv->get_fn_order() +
          fw->get_fn_order() + ru->get_inv_ref_order();

  limit_order(o);
  fu->set_quad_order(o, FN_VAL);
  fv->set_quad_order(o);
  fw->set_quad_order(o, FN_VAL);

  double* uval = fu->get_fn_values();
  double* wval = fw->get_fn_values();
  double *dvdx, *dvdy;
  fv->get_dx_dy_values(dvdx, dvdy); // 'v' solution => derivatives already transformed to physical element

  double result = 0.0;
  h1_integrate_dd_expression(uval[i] * dvdx[i] * wval[i]);
  return result;
}

inline double int_u_dvdy_w(RealFunction* fu, RealFunction* fv, RealFunction* fw, RefMap* ru, RefMap* rv, RefMap* rw)
// assumes that 'v' is a solution
{
  Quad2D* quad = fu->get_quad_2d();

  int o = fu->get_fn_order() + fv->get_fn_order() +
          fw->get_fn_order() + ru->get_inv_ref_order();

  limit_order(o);
  fu->set_quad_order(o, FN_VAL);
  fv->set_quad_order(o);
  fw->set_quad_order(o, FN_VAL);

  double* uval = fu->get_fn_values();
  double* wval = fw->get_fn_values();
  double *dvdx, *dvdy;
  fv->get_dx_dy_values(dvdx, dvdy); // 'v' solution => derivatives already transformed to physical element

  double result = 0.0;
  h1_integrate_dd_expression(uval[i] * dvdy[i] * wval[i]);
  return result;
}


// definition of boundary conditions
int xvel_bc_type(int marker) {
  if (marker == 2) return BC_NONE;
  else return BC_ESSENTIAL;
}

scalar xvel_bc_value(int marker, double x, double y) {
  if (marker == 4) {
    // time-dependent inlet velocity
    //double val_y = VEL_INLET; //constant profile
    double val_y = VEL_INLET * y*(H-y) / (H/2.)/(H/2.); //parabolic profile with peak VEL_INLET at y = H/2
    if (TIME <= STARTUP_TIME) return val_y * TIME/STARTUP_TIME;
    else return val_y;
  }
  else return 0;
}

int yvel_bc_type(int marker) {
  if (marker == 2) return BC_NONE;
  else return BC_ESSENTIAL;
}

int press_bc_type(int marker)
  { return BC_NONE; }

// velocities from the previous time step and for the Newton's iteration
Solution xprev, yprev, pprev, xiter, yiter, piter;

// bilinear and linear forms corresponding to simple linearization
// of convective term
scalar bilinear_form_sym_0_0_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv) / RE +
           int_u_v(fu, fv, ru, rv) / TAU; }

scalar simple_bilinear_form_unsym_0_0_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv); }

scalar simple_linear_form_0(RealFunction* fv, RefMap* rv)
  { return int_u_v(&xprev, fv, xprev.get_refmap(), rv) / TAU; }

scalar simple_linear_form_1(RealFunction* fv, RefMap* rv)
  { return int_u_v(&yprev, fv, yprev.get_refmap(), rv) / TAU; }

scalar bilinear_form_unsym_0_2(RealFunction* fp, RealFunction* fv, RefMap* rp, RefMap* rv)
  { return -int_u_dvdx(fp, fv, rp, rv); }

scalar bilinear_form_unsym_1_2(RealFunction* fp, RealFunction* fv, RefMap* rp, RefMap* rv)
  { return -int_u_dvdy(fp, fv, rp, rv); }

// bilinear and linear forms corresponding to the Newton's method
scalar newton_bilinear_form_unsym_0_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return   int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv)
           + int_u_dvdx_w(fu, &xprev, fv, ru, xprev.get_refmap(), rv); }

scalar newton_bilinear_form_unsym_0_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_u_dvdy_w(fu, &xprev, fv, ru, xprev.get_refmap(), rv); }

scalar newton_bilinear_form_unsym_1_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_u_dvdx_w(fu, &yprev, fv, ru, yprev.get_refmap(), rv); }

scalar newton_bilinear_form_unsym_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return   int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv)
           + int_u_dvdy_w(fu, &yprev, fv, ru, yprev.get_refmap(), rv); }

scalar newton_F_0(RealFunction* fv, RefMap* rv)
  { return   int_u_v(&xiter, fv, xiter.get_refmap(), rv) / TAU
           - int_u_v(&xprev, fv, xprev.get_refmap(), rv) / TAU
           + int_grad_u_grad_v_II(&xiter, fv, xiter.get_refmap(), rv) / RE   // The 'II' version of int_grad_u_grad_v()
             // does exactly the same as the original function, but it assumes that the first argument is a Solution.
             // This inconsistency has historical reasons and we are going to eliminate it soon.
           + int_w_nabla_u_v_II(&xiter, &yiter, &xiter, fv, xiter.get_refmap(), rv)
           - int_u_dvdx(&piter, fv, piter.get_refmap(), rv); }

scalar newton_F_1(RealFunction* fv, RefMap* rv)
  { return   int_u_v(&yiter, fv, yiter.get_refmap(), rv) / TAU
           - int_u_v(&yprev, fv, yprev.get_refmap(), rv) / TAU
           + int_grad_u_grad_v_II(&yiter, fv, yiter.get_refmap(), rv) / RE
           + int_w_nabla_u_v_II(&xiter, &yiter, &yiter, fv, yiter.get_refmap(), rv)
           - int_u_dvdy(&piter, fv, piter.get_refmap(), rv); }

scalar newton_F_2(RealFunction* fv, RefMap* rv)
  { return   int_dudx_v_II(&xiter, fv, xiter.get_refmap(), rv)
           + int_dudy_v_II(&yiter, fv, yiter.get_refmap(), rv); }


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");

  // a-priori mesh refinements
  mesh.refine_element(3, 2);
  mesh.refine_element(0, 0);
  mesh.refine_element(1, 1);
  mesh.refine_element(2, 0);
  mesh.refine_element(4, 2);
  mesh.refine_element(5, 2);
  mesh.refine_element(6, 2);
  mesh.refine_element(7, 2);
  mesh.refine_element(8, 2);
  mesh.refine_element(9, 0);
  mesh.refine_element(10, 1);
  mesh.refine_element(11, 0);
  mesh.refine_element(12, 2);
  mesh.refine_element(13, 2);
  mesh.refine_element(14, 0);
  mesh.refine_element(15, 0);
  mesh.refine_element(26, 0);
  mesh.refine_element(27, 0);
  mesh.refine_element(32, 2);
  mesh.refine_element(33, 2);
  mesh.refine_element(34, 2);
  mesh.refine_element(35, 2);
  mesh.refine_element(46, 0);
  mesh.refine_element(47, 0);
  mesh.refine_element(48, 0);
  mesh.refine_element(49, 0);
  mesh.refine_all_elements();
  mesh.refine_towards_boundary(5, 4, false);
  mesh.refine_towards_boundary(1, 4);
  mesh.refine_towards_boundary(3, 4);

  // display the mesh
  //MeshView mview("Hello world!", 100, 100, 1100, 400);
  //mview.show(&mesh);
  //mview.wait_for_keypress();

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

  // initial condition: xprev and yprev are zero
  xprev.set_zero(&mesh);
  yprev.set_zero(&mesh);
  pprev.set_zero(&mesh);

  // set up weak formulation
  WeakForm wf(3);
  if (NEWTON_ON) {
    wf.add_biform(0, 0, bilinear_form_sym_0_0_1_1, SYM);
    wf.add_biform(0, 0, newton_bilinear_form_unsym_0_0, UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(0, 1, newton_bilinear_form_unsym_0_1, UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(0, 2, bilinear_form_unsym_0_2, ANTISYM);
    wf.add_biform(1, 0, newton_bilinear_form_unsym_1_0, UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(1, 1, bilinear_form_sym_0_0_1_1, SYM);
    wf.add_biform(1, 1, newton_bilinear_form_unsym_1_1, UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(1, 2, bilinear_form_unsym_1_2, ANTISYM);
    wf.add_liform(0, newton_F_0, ANY, 5, &xprev, &yprev, &xiter, &yiter, &piter);
    wf.add_liform(1, newton_F_1, ANY, 5, &xprev, &yprev, &xiter, &yiter, &piter);
    wf.add_liform(2, newton_F_2, ANY, 2, &xiter, &yiter);
  }
  else {
    wf.add_biform(0, 0, bilinear_form_sym_0_0_1_1, SYM);
    wf.add_biform(0, 0, simple_bilinear_form_unsym_0_0_1_1, UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(1, 1, bilinear_form_sym_0_0_1_1, SYM);
    wf.add_biform(1, 1, simple_bilinear_form_unsym_0_0_1_1, UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(0, 2, bilinear_form_unsym_0_2, ANTISYM);
    wf.add_biform(1, 2, bilinear_form_unsym_1_2, ANTISYM);
    wf.add_liform(0, simple_linear_form_0, ANY, 1, &xprev);
    wf.add_liform(1, simple_linear_form_1, ANY, 1, &yprev);
  }

  // visualization
  VectorView vview("velocity [m/s]", 0, 0, 1500, 470);
  ScalarView pview("pressure [Pa]", 0, 530, 1500, 470);
  vview.set_min_max_range(0, 1.6);
  //pview.set_min_max_range(-0.9, 1.0);
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
      // initialize the Newton's iteration
      xiter.copy(&xprev);
      yiter.copy(&yprev);
      piter.copy(&pprev);

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
      pprev.copy(&piter);
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

  // wait for keypress or mouse input
  View::wait();
}
