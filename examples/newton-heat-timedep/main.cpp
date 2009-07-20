#include "hermes2d.h"
#include "solver_umfpack.h"

//
//  This example demonstrates the employment of the Newton's method to
//  a nonlinear time-dependent PDE discretized implicitly in time
//  (via implicit Euler or Crank-Nicolson). Some problem parameters can
//  be changed below.
//
//  PDE: non-stationary heat transfer with nonlinear thermal conductivity
//  HEATCAP*dT/dt - div[lambda(T)grad T] = 0
//
//  Domain: square
//
//  BC:  T = 100 on the left, top and bottom edges
//       dT/dn = 0 on the right edge
//
//  Time-stepping: either implicit Euler or Crank-Nicolson
//

/********** Problem parameters ***********/

int TIME_DISCR = 2;        // 1 for implicit Euler, 2 for Crank-Nicolson
int PROJ_TYPE = 1;         // 1 for H1 projections, 0 for L2 projections
double HEATCAP = 1e6;      // heat capacity
double TAU = 0.5;          // time step
int NSTEP = 1000;          // number of time steps to do
double NEWTON_TOL = 1e-3;  // convergence criterion for the Newton's method

// thermal conductivity (temperature-dependent
// for any u, this function has to be  positive in the entire domain!
double lam(double T)  { return 10 + 0.1*pow(T, 2); }
double dlam_dT(double T) { return 0.1*2*pow(T, 1); }

/********** Definition of boundary conditions ***********/

int bc_types(int marker)
{
 if (marker == 4 || marker == 1 || marker == 3) return BC_ESSENTIAL;
//   if (marker == 4) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

scalar bc_values(int marker, double x, double y)
{
 if (marker == 4 || marker == 1 || marker == 3) return 100;
//   if (marker == 4) return -4.0 * sqr(y) + 4.0 * y;
  else return 0.0;
}


/********** Definition of Jacobian matrices and residual vectors ***********/

// Residuum for the Euler time discretization
inline double F_euler(RealFunction* Tprev, RealFunction* Titer, RealFunction* fu, RefMap* ru)
{
  Quad2D* quad = fu->get_quad_2d();
  RefMap* rv = ru;

  int o = 3 * Titer->get_fn_order() + fu->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Tprev->set_quad_order(o, FN_VAL);
  Titer->set_quad_order(o);
  fu->set_quad_order(o);

  double* Titer_val = Titer->get_fn_values();
  double* Tprev_val = Tprev->get_fn_values();
  double* uval = fu->get_fn_values();
  double *dTiter_dx, *dTiter_dy, *dudx, *dudy;
  Titer->get_dx_dy_values(dTiter_dx, dTiter_dy);
  fu->get_dx_dy_values(dudx, dudy);

  // u is a test function
  double result;
  h1_integrate_dd_expression(( HEATCAP*(Titer_val[i] - Tprev_val[i])*uval[i]/TAU +
                               lam(Titer_val[i]) * (dTiter_dx[i]*t_dudx + dTiter_dy[i]*t_dudy)));

  return result;
}

// Jacobian matrix for the implicit Euler time discretization
inline double J_euler(RealFunction* Titer, RealFunction* fu,
                RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 2 * Titer->get_fn_order() + fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Titer->set_quad_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double* Titer_val = Titer->get_fn_values();
  double* uval = fu->get_fn_values();
  double* vval = fv->get_fn_values();

  double *dTiter_dx, *dTiter_dy, *dudx, *dudy, *dvdx, *dvdy;
  Titer->get_dx_dy_values(dTiter_dx, dTiter_dy);
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);

  // u is a basis function, v a test function
  double result;
  h1_integrate_dd_expression(( HEATCAP * uval[i] * vval[i] / TAU +
                               dlam_dT(Titer_val[i]) * uval[i] * (dTiter_dx[i]*t_dvdx + dTiter_dy[i]*t_dvdy) +
                               lam(Titer_val[i]) * (t_dudx*t_dvdx + t_dudy*t_dvdy)));

  return result;
}

// Residuum for the Crank-Nicolson time discretization
inline double F_cranic(RealFunction* Tprev, RealFunction* Titer, RealFunction* fu, RefMap* ru)
{
  Quad2D* quad = fu->get_quad_2d();
  RefMap* rv = ru;

  int o = 3 * Titer->get_fn_order() + fu->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Tprev->set_quad_order(o);
  Titer->set_quad_order(o);
  fu->set_quad_order(o);

  double* Titer_val = Titer->get_fn_values();
  double* Tprev_val = Tprev->get_fn_values();
  double* uval = fu->get_fn_values();
  double *dTiter_dx, *dTiter_dy, *dTprev_dx, *dTprev_dy, *dudx, *dudy;
  Titer->get_dx_dy_values(dTiter_dx, dTiter_dy);
  Tprev->get_dx_dy_values(dTprev_dx, dTprev_dy);
  fu->get_dx_dy_values(dudx, dudy);

  // u is a test function
  double result;
  h1_integrate_dd_expression(( HEATCAP * (Titer_val[i] - Tprev_val[i]) * uval[i] / TAU +
                               0.5 * lam(Titer_val[i]) * (dTiter_dx[i]*t_dudx + dTiter_dy[i]*t_dudy) +
                               0.5 * lam(Tprev_val[i]) * (dTprev_dx[i]*t_dudx + dTprev_dy[i]*t_dudy)
                            ));

  return result;
}

// Jacobian matrix for the Crank-Nicolson time discretization
inline double J_cranic(RealFunction* Titer, RealFunction* fu,
                RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 2 * Titer->get_fn_order() + fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Titer->set_quad_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double* Titer_val = Titer->get_fn_values();
  double* uval = fu->get_fn_values();
  double* vval = fv->get_fn_values();

  double *dTiter_dx, *dTiter_dy, *dudx, *dudy, *dvdx, *dvdy;
  Titer->get_dx_dy_values(dTiter_dx, dTiter_dy);
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);

  // u is a basis function, v a test function
  double result;
  h1_integrate_dd_expression(( HEATCAP * uval[i] * vval[i] / TAU +
                               0.5 * dlam_dT(Titer_val[i]) * uval[i] * (dTiter_dx[i]*t_dvdx + dTiter_dy[i]*t_dvdy) +
                               0.5 * lam(Titer_val[i]) * (t_dudx*t_dvdx + t_dudy*t_dvdy)
                            ));

  return result;
}

/********** Definition of linear and bilinear forms for Hermes ***********/

Solution Tprev, // previous time step solution, for the time integration method
         Titer; // solution converging during the Newton's iteration

// Implicit Euler method (1st-order in time)
scalar linear_form_0_euler(RealFunction* fv, RefMap* rv)
{ return F_euler(&Tprev, &Titer, fv, rv); }
scalar bilinear_form_0_0_euler(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ return J_euler(&Titer, fu, fv, ru, rv); }

// Crank-Nicolson method (2nd-order in time)
scalar linear_form_0_cranic(RealFunction* fv, RefMap* rv)
{ return F_cranic(&Tprev, &Titer, fv, rv); }
scalar bilinear_form_0_0_cranic(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ return J_cranic(&Titer, fu, fv, ru, rv); }

// *************************************************************

int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("square.mesh");
  for(int i = 0; i < 5; i++) mesh.refine_all_elements();

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(1);
  space.assign_dofs();

  WeakForm wf(1);
  if(TIME_DISCR == 1) {
    wf.add_biform(0, 0, bilinear_form_0_0_euler, UNSYM, ANY, 1, &Titer);
    wf.add_liform(0, linear_form_0_euler, ANY, 2, &Titer, &Tprev);
  }
  else {
    wf.add_biform(0, 0, bilinear_form_0_0_cranic, UNSYM, ANY, 1, &Titer);
    wf.add_liform(0, linear_form_0_cranic, ANY, 2, &Titer, &Tprev);
  }

  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  char title[100];
  ScalarView view("", 0, 0, 600, 600);
  ScalarView view2("", 700, 0, 600, 600);

  // setting the Dirichlet lift to be the initial condition
  Tprev.set_dirichlet_lift(&space, &pss);
  nls.set_ic(&Tprev, &Tprev, PROJ_TYPE);
  Titer.copy(&Tprev);

  // view initial guess for Newton's method
  sprintf(title, "Initial guess for the Newton's method");
  view.set_title(title);
  view.show(&Titer);
  printf("Click into the image window and press any key.\n");
  view.wait_for_keypress();

  Solution sln;
  // time stepping
  for(int n = 1; n <= NSTEP; n++)
  {

    info("\n---- Time step %d -----------------------------------------------", n);

    // set initial condition for the Newton's iteration
    // actually needed only when space changes
    // otherwise initial solution vector is that one
    // from the previous time level
    //nls.set_ic(&Titer, &Titer);

    int it = 1;
    double res_l2_norm;
    do
    {
      info("\n---- Time step %d, Newton iter %d ---------------------------------\n", n, it++);

      nls.assemble();
      nls.solve(1, &sln);

      res_l2_norm = nls.get_residuum_l2_norm();
      info("Residuum L2 norm: %g\n", res_l2_norm);
      // want to see Newtons iterations
//       sprintf(title, "Time level %d, Newton iteration %d", n, it-1);
//       view.set_title(title);
//       view.show(&sln);
//       view.wait_for_keypress();

      Titer = sln;

    }
    while (res_l2_norm > NEWTON_TOL);

    // visualization of solution on the n-th time level
    sprintf(title, "Time level %d", n);
    //view.set_min_max_range(90,100);
    view.set_title(title);
    view.show(&Titer);
    //view.wait_for_keypress();

    // uncomment one of the following lines to generate a series of video frames
    //vview.save_numbered_screenshot("velocity%03d.bmp", i, true);
    //pview.save_numbered_screenshot("pressure%03d.bmp", i, true);
    // the frames can then be converted to a video file with the command
    // mencoder "mf://velocity*.bmp" -mf fps=20 -o velocity.avi -ovc lavc -lavcopts vcodec=mpeg4



    // copying result of the Newton's iteration into Tprev
    Tprev.copy(&Titer);
  }

  printf("Click into the image window and press 'q' to finish.\n");
  View::wait();
  return 0;
}
