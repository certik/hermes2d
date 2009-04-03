#include "hermes2d.h"
#include "solver_umfpack.h"

// 
//  Nonlinear solver test:
//
//  PDE: non stationary heat transfer with nonlinear thermal conductuvity
//  dT/dt - div[lambda(T)grad T] = 0
//
//  Domain: square
//
//  BC:  T = 100 on the left, top and bottom edges
//       dT/dn = 0 on the right edge
//

// time step and number of time steps
double TAU = 1e-7;
int NSTEP = 200;

// thermal conductivity (temperature-dependent)
// for any u, this function has to be  positive in the entire domain!
double lam(double T) { return 1.0 + T*T; } 
double dlam_dT(double T) { return 2*T; }


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

Solution Tprev, // previous time step solution, for the implicit Euler method
         Titer; // solution converging during the Newton's iteration

inline double F(RealFunction* Tprev, RealFunction* Titer, RealFunction* fu, RefMap* ru)
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

  h1_integrate_dd_expression(( (Titer_val[i] - Tprev_val[i])*uval[i]/TAU + 
                               lam(Titer_val[i]) * (dTiter_dx[i]*t_dudx + dTiter_dy[i]*t_dudy)));

  return result;
}

scalar linear_form_0(RealFunction* fv, RefMap* rv)
{ return F(&Tprev, &Titer, fv, rv); }

inline double J(RealFunction* Titer, RealFunction* fu, 
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

  h1_integrate_dd_expression(( vval[i] * uval[i]/TAU +
                               dlam_dT(Titer_val[i]) * uval[i] * (dTiter_dx[i]*t_dvdx + dTiter_dy[i]*t_dvdy) +
                               lam(Titer_val[i]) * (t_dudx*t_dvdx + t_dudy*t_dvdy)));

  return result;
}

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
  wf.add_liform(0, linear_form_0, ANY, 2, &Titer, &Tprev);
  
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  char title[100];  
  ScalarView view("", 0, 0, 600, 600);
  ScalarView view2("", 700, 0, 600, 600);
  
  // initial condition at zero time level
  Tprev.set_const(&mesh, 0.0);
  nls.set_ic(&Tprev, &Tprev);
  Titer.copy(&Tprev);
//   Tprev.set_const(&mesh, 0.0);
//   Titer.copy(&Tprev);
//   nls.set_ic(&Titer, &Titer);

  // view initial guess for Newton's method
  // satisfies BC conditions 
  sprintf(title, "Initial iteration");
  view.set_title(title);
  view.show(&Titer);    
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
    nls.set_ic(&Titer, &Titer);
 
    int it = 1; 
    double res_l2_norm; 
    do
    {
      info("\n---- Time step %d, iter %d ---------------------------------------\n", n, it++);
    
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
    while (res_l2_norm > 1e-4);

    // visualization of solution on the n-th time level
    sprintf(title, "Time level %d", n);
    view.set_min_max_range(90,100);
    view.set_title(title);
    view.show(&Titer);    
    //view.wait_for_keypress();

    // copying result of the Newton's iteration into Tprev
    Tprev.copy(&Titer);
  }  

  View::wait();
  return 0;
}
