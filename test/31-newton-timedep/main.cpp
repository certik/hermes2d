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
  
  double result = 0.0;
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

  double result = 0.0;
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
  Mesh mesh, basemesh;
  basemesh.load("square.mesh");
  for(int i = 0; i < 2; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);
  
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(2);
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
  OrderView ordview("Polynomial Orders", 0, 700, 600, 600);
  
  GnuplotGraph graph_err;
  graph_err.set_captions("","Time step","Error");
  graph_err.add_row();
  GnuplotGraph graph_dofs;
  graph_dofs.set_captions("","Time step","DOFs");
  graph_dofs.add_row();
  
  // initial condition at zero time level
  Tprev.set_const(&mesh, 0.0);
  nls.set_ic(&Tprev, &Tprev);
  Titer.copy(&Tprev);

  // view initial guess for Newton's method
  // satisfies BC conditions 
  sprintf(title, "Initial iteration");
  view.set_title(title);
  view.show(&Titer);    
  view.wait_for_keypress();

  Solution sln, rsln;
 // time stepping
  for(int n = 1; n <= NSTEP; n++) 
  {

    info("\n---- Time step %d -----------------------------------------------", n);
    
    // unrefinements
    if (n % 10 == 2) {
      mesh.copy(&basemesh);
      space.set_uniform_order(2);
    }
    
    int at = 0;
    bool done = false;
    double err;
    do 
    {
      info("\n---- Time step %d, Adaptivity %d ---------------------------------------\n", n, ++at);

    
      int it = 1; 
      double res_l2_norm; 
      space.assign_dofs();    
      if (n > 1 || at > 1) nls.set_ic(&rsln, &Titer);
      else nls.set_ic(&Titer, &Titer);
      do
      {
        info("\n---- Time step %d, Adaptivity %d, iter %d ---------------------------------------\n", n, at, it++);
        
        nls.assemble();
        nls.solve(1, &sln);
  
        res_l2_norm = nls.get_residuum_l2_norm(); 
        info("Residuum L2 norm: %g\n", res_l2_norm);
   
        Titer.copy(&sln);          
      }
      while (res_l2_norm > 1e-4);
      
      ordview.show(&space);
      
      // reference solution
      it = 1;
      RefNonlinSystem rs(&nls);
      rs.prepare();
      if (n > 1 || at > 1) rs.set_ic(&rsln, &Titer);
      else rs.set_ic(&Titer, &Titer);
      do
      {
        info("\n---- Time step %d, Adaptivity %d, iter %d (Reference solution) ---------------------------------------\n", n, at, it++);
        
        rs.assemble();
        rs.solve(1, &rsln);

        res_l2_norm = rs.get_residuum_l2_norm(); 
        info("Residuum L2 norm: %g\n", res_l2_norm);
  
        Titer.copy(&rsln);          
      }
      while (res_l2_norm > 1e-4);

      H1OrthoHP hp(1, &space);
      err = hp.calc_error(&sln, &rsln) * 100;
      info("Error: %g", err);
      if (err < 0.001) done = true;
      else hp.adapt(0.3, 1);

    }
    while (!done);
      
    // visualization of solution on the n-th time level
    sprintf(title, "Time level %d", n);
    view.set_min_max_range(90,100);
    view.set_title(title);
    view.show(&Titer);    
    //view.wait_for_keypress();

    graph_err.add_values(0, n, err);
    graph_err.save("error.txt");
    graph_dofs.add_values(0, n, space.get_num_dofs());
    graph_dofs.save("dofs.txt");
    
    // copying result of the Newton's iteration into Tprev
    Tprev.copy(&Titer);
  }  

  View::wait();
  return 0;
}
