#include "hermes2d.h"
#include "solver_umfpack.h"

// 
//  PDE: non-stationary complex Schroedinger equation 
//  ih \partial \psi/\partial t = -h^2/(2m) \Delta \psi +
//  g \psi |\psi|^2 + 1/2 m \omega^2 (x^2 + y^2) \psi
//
//  Domain: square of edge length 2 centered at origin
//
//  BC:  homogeneous Dirichlet everywhere on th eboundary
//
//  Time-stepping: either implicit Euler or Crank-Nicolson
//

//#include <complex.h>

/********** Problem parameters ***********/ 
double H = 1;         //Planck constant 6.626068e-34;
double M = 1; 
double G = 1; 
double OMEGA = 1;     
int TIME_DISCR = 1;   // 1 for implicit Euler, 2 for Crank-Nicolson
int PROJ_TYPE = 2;    // 1 for H1 projections, 2 for L2 projections
double TAU = 0.5;     // time step
int NSTEP = 100;      // number of time steps to do

/********** Definition of initial conditions ***********/ 

scalar fn_init(double x, double y, scalar& dx, scalar& dy) {
  return exp(-x*x - y*y);
}

/********** Definition of boundary conditions ***********/ 

int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

complex bc_values(int marker, double x, double y)
{
  return complex(0.0, 0.0);
}

/********** Definition of Jacobian matrices and residual vectors ***********/ 

// Residuum for the Euler time discretization
inline complex F_euler(ScalarFunction* Psi_prev, ScalarFunction* Psi_iter, RealFunction* fu, RefMap* ru)
{
  scalar ii = complex(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Quad2D* quad = fu->get_quad_2d();
  RefMap* rv = ru;

  //not sure here:
  int o = 3 * Psi_iter->get_fn_order() + fu->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Psi_prev->set_quad_order(o, FN_VAL);
  Psi_iter->set_quad_order(o);
  fu->set_quad_order(o);

  scalar* Psi_iter_val = Psi_iter->get_fn_values();
  scalar* Psi_prev_val = Psi_prev->get_fn_values();
  scalar* dPsi_iter_dx, *dPsi_iter_dy;
  double* uval = fu->get_fn_values();
  double *dudx, *dudy;
  Psi_iter->get_dx_dy_values(dPsi_iter_dx, dPsi_iter_dy);
  fu->get_dx_dy_values(dudx, dudy);

  // obtain physical coordinates of int. points
  double* x = ru->get_phys_x(o);  
  double* y = ru->get_phys_y(o);  

  // u is a test function
  h1_integrate_dd_expression(( 
    ii * H * (Psi_iter_val[i] - Psi_prev_val[i]) * uval[i] / TAU 
    - H*H/(2*M) * (dPsi_iter_dx[i]*t_dudx + dPsi_iter_dy[i]*t_dudy)
    - G * Psi_iter_val[i] *  Psi_iter_val[i] * conj(Psi_iter_val[i]) * uval[i]
    - .5*M*OMEGA*OMEGA * (x[i]*x[i] + y[i]*y[i]) * Psi_iter_val[i] * uval[i]
  ));

  return result;
}

// Jacobian matrix for the implicit Euler time discretization
inline complex J_euler(ScalarFunction* Psi_iter, RealFunction* fu, 
                RealFunction* fv, RefMap* ru, RefMap* rv)
{
  scalar ii = complex(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Quad2D* quad = fu->get_quad_2d();

  int o = 2 * Psi_iter->get_fn_order() + fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Psi_iter->set_quad_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  scalar* Psi_iter_val = Psi_iter->get_fn_values();
  scalar *dPsi_iter_dx, *dPsi_iter_dy;
  double* uval = fu->get_fn_values();
  double* vval = fv->get_fn_values();
  double *dudx, *dudy, *dvdx, *dvdy;
  Psi_iter->get_dx_dy_values(dPsi_iter_dx, dPsi_iter_dy);
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);

  // obtain physical coordinates of int. points
  double* x = ru->get_phys_x(o);  
  double* y = ru->get_phys_y(o);  

  // u is a basis function, v a test function
  h1_integrate_dd_expression(( 
    ii * H * uval[i] * vval[i] / TAU 
    - H*H/(2*M) * (t_dudx*t_dvdx + t_dudy*t_dvdy)
    - 2* G * uval[i] *  Psi_iter_val[i] * conj(Psi_iter_val[i]) * vval[i]
    - G * Psi_iter_val[i] *  Psi_iter_val[i] * conj(uval[i]) * vval[i]
    - .5*M*OMEGA*OMEGA * (x[i]*x[i] + y[i]*y[i]) * uval[i] * vval[i]
  ));

  /* old code from the heat transfer equation
  h1_integrate_dd_expression(( HEATCAP * uval[i] * vval[i] / TAU +
                               dlam_dT(Psi_iter_val[i]) * uval[i] * (dPsi_iter_dx[i]*t_dvdx + 
                               dPsi_iter_dy[i]*t_dvdy) +
                               lam(Psi_iter_val[i]) * (t_dudx*t_dvdx + t_dudy*t_dvdy)));
  */
  return result;
}

/********** Definition of linear and bilinear forms for Hermes ***********/ 

Solution Psi_prev, // previous time step solution, for the time integration method
         Psi_iter; // solution converging during the Newton's iteration

// Implicit Euler method (1st-order in time)
complex linear_form_0_euler(RealFunction* fv, RefMap* rv)
{ return F_euler(&Psi_prev, &Psi_iter, fv, rv); }
complex bilinear_form_0_0_euler(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ return J_euler(&Psi_iter, fu, fv, ru, rv); }

/*
// Crank-Nicolson method (2nd-order in time)
scalar linear_form_0_cranic(RealFunction* fv, RefMap* rv)
{ return F_cranic(&Psi_prev, &Psi_iter, fv, rv); }
scalar bilinear_form_0_0_cranic(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ return J_cranic(&Psi_iter, fu, fv, ru, rv); }
*/

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
    wf.add_biform(0, 0, bilinear_form_0_0_euler, UNSYM, ANY, 1, &Psi_iter);
    wf.add_liform(0, linear_form_0_euler, ANY, 2, &Psi_iter, &Psi_prev);
  }
  else {
    printf("Crank-Nicolson not implemented yet.\n");
    exit(0);
    /*
    wf.add_biform(0, 0, bilinear_form_0_0_cranic, UNSYM, ANY, 1, &Psi_iter);
    wf.add_liform(0, linear_form_0_cranic, ANY, 2, &Psi_iter, &Psi_prev);
    */
  }

  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  char title[100];  
  ScalarView view("", 0, 0, 600, 600);
  ScalarView view2("", 700, 0, 600, 600);
  
  // setting initial condition at zero time level
  Psi_prev.set_exact(&mesh, fn_init);
  nls.set_ic(&Psi_prev, &Psi_prev, PROJ_TYPE);
  Psi_iter.copy(&Psi_prev);

  // view initial guess for Newton's method 
  //sprintf(title, "Initial guess for the Newton's method");
  //view.set_title(title);
  //view.show(&Psi_iter);    
  //view.wait_for_keypress();

  Solution sln;
  // time stepping
  for(int n = 1; n <= NSTEP; n++) 
  {

    info("\n---- Time step %d -----------------------------------------------", n);
    
    // set initial condition for the Newton's iteration
    // actually needed only when space changes
    // otherwise initial solution vector is that one 
    // from the previous time level
    //nls.set_ic(&Psi_iter, &Psi_iter);
 
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

      Psi_iter = sln;
        
    }
    while (res_l2_norm > 1e-4);

    // visualization of solution on the n-th time level
    sprintf(title, "Time level %d", n);
    //view.set_min_max_range(90,100);
    view.set_title(title);
    view.show(&Psi_iter);    
    //view.wait_for_keypress();

    // copying result of the Newton's iteration into Psi_prev
    Psi_prev.copy(&Psi_iter);
  }  

  View::wait();
  return 0;
}
