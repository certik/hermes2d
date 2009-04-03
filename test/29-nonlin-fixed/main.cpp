// 
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

//#define DEBUG_ORDER

#include "hermes2d.h"
#include "solver_umfpack.h"

// time step and number of time steps
double TAU = 1e-7;
int NSTEP = 50;

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

inline double int_lambda_grad_u_grad_v(RealFunction* w, RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order() + 2*w->get_fn_order();
  limit_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);
  w->set_quad_order(o);

  double *dudx, *dudy, *dvdx, *dvdy;
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);
  double* wval = w->get_fn_values();
  
  h1_integrate_dd_expression(lam(wval[i]) * (t_dudx * t_dvdx + t_dudy * t_dvdy));
  return result;
}


Solution Tprev, Titer;

scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_lambda_grad_u_grad_v(&Titer, fu, fv, ru, rv) + int_u_v(fu, fv, ru, rv) / TAU;
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return int_w_v(&Tprev, fv, rv) / TAU;
}
 

int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("square.mesh");
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(1);
  space.assign_dofs();
  
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form, UNSYM, ANY, 1, &Titer);
  wf.add_liform(0, linear_form, ANY, 1, &Tprev);
  
  UmfpackSolver umfpack;
  LinSystem ls(&wf, &umfpack);
  ls.set_spaces(1, &space);
  ls.set_pss(1, &pss);
  
  ScalarView view("Iteration", 0, 0, 880, 800);
  
  Tprev.set_zero(&mesh);
  Titer.set_zero(&mesh);
  

  for(int n = 1; n <= NSTEP; n++) 
  {

    info("\n---- Time step %d -----------------------------------------------", n);

    double residuum, error;  
    int it = 1;
    do
    {
      info("------------------ it = %d------------------------------", it++);
  
      Solution sln;
      ls.assemble();  
      ls.solve(1, &sln);
      info("residuum: %g", residuum = l2_error(&Titer, &sln)); 
      
  
      view.show(&sln);
      view.wait_for_keypress();
  
      Titer = sln;
    }
    while (residuum > 1e-4);
    
    Tprev.copy(&Titer);
  }
  
  
  View::wait();
  return 0;
}
