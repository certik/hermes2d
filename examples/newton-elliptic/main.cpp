#include "hermes2d.h"
#include "solver_umfpack.h"
#include "function.h"

//
//  The purpose of this example is to show how to use the Newton's
//  method for a stationary nonlinear PDE problem. Some problem
//  parameters can be changed below.
//
//  PDE: stationary heat transfer with nonlinear thermal conductivity
//  - div[lambda(T)grad T] = 0
//
//  Domain: square
//
//  BC:  T = 100 on the left, top and bottom edges
//       dT/dn = 0 on the right edge
//

/********** Problem parameters ***********/

int PROJ_TYPE = 0;         // 1 for H1 projections, 0 for L2 projections
double NEWTON_TOL = 1e-3;  // convergence criterion for the Newton's method

// thermal conductivity (temperature-dependent)
// for any u, this function has to be  positive in the entire domain!
double lam(double T) { return 10 + 0.1*pow(T, 2); }
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

// Residuum
inline double F(RealFunction* Titer, RealFunction* fu, RefMap* ru)
{
  Quad2D* quad = fu->get_quad_2d();
  RefMap* rv = ru;

  int o = 3 * Titer->get_fn_order() + fu->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Titer->set_quad_order(o);
  fu->set_quad_order(o);

  double* Titer_val = Titer->get_fn_values();
  double* uval = fu->get_fn_values();
  double *dTiter_dx, *dTiter_dy, *dudx, *dudy;
  Titer->get_dx_dy_values(dTiter_dx, dTiter_dy);
  fu->get_dx_dy_values(dudx, dudy);

  // u is a test function
  double result;
  h1_integrate_dd_expression(( lam(Titer_val[i]) * (dTiter_dx[i]*t_dudx + dTiter_dy[i]*t_dudy)));

  return result;
}

// Jacobian matrix
inline double J(RealFunction* Titer, RealFunction* fu,
                RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 2 * Titer->get_fn_order() + fu->get_fn_order()
     + fv->get_fn_order() + ru->get_inv_ref_order();
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
  h1_integrate_dd_expression(( dlam_dT(Titer_val[i]) * uval[i]
    * (dTiter_dx[i]*t_dvdx + dTiter_dy[i]*t_dvdy) +
    lam(Titer_val[i]) * (t_dudx*t_dvdx + t_dudy*t_dvdy)));

  return result;
}

/********** Definition of linear and bilinear forms for Hermes ***********/

Solution Titer; // solution converging during the Newton's iteration

scalar linear_form_0(RealFunction* fv, RefMap* rv)
{ return F(&Titer, fv, rv); }
scalar bilinear_form_0_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ return J(&Titer, fu, fv, ru, rv); }

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
  wf.add_biform(0, 0, bilinear_form_0_0, UNSYM, ANY, 1, &Titer);
  wf.add_liform(0, linear_form_0, ANY, 1, &Titer);

  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  char title[100];
  ScalarView view("", 0, 0, 600, 600);
  ScalarView view2("", 700, 0, 600, 600);

  // setting the Dirichlet lift to be the initial condition
  Titer.set_dirichlet_lift(&space, &pss);
  nls.set_ic(&Titer, &Titer, PROJ_TYPE);

  // view initial guess for Newton's method
  /*
  sprintf(title, "Initial guess for the Newton's method");
  view.set_title(title);
  view.show(&Titer);
  view.wait_for_keypress();
  */
  Solution sln;

  int it = 1;
  double res_l2_norm;
  do
  {
    info("\n---- Newton iter %d ---------------------------------\n", it++);

    nls.assemble();
    nls.solve(1, &sln);

    res_l2_norm = nls.get_residuum_l2_norm();
    info("Residuum L2 norm: %g\n", res_l2_norm);
    // want to see Newtons iterations
    sprintf(title, "Newton iteration %d", it-1);
    view.set_title(title);
    view.show(&sln);
    printf("Click into the image window and press any key to proceed.\n");
    view.wait_for_keypress();

    Titer = sln;

  }
  while (res_l2_norm > NEWTON_TOL);

  //View::wait();
  return 0;
}

