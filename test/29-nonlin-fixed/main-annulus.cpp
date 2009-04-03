// 
//  Nonlinear solver test:
//
//  PDE:  -\Delta u + u^3 = 0
//
//  Domain: annulus centered at origin, with inner radius A = 0.05 and
//  outer radius B = 0.5.
//
//  BC:  u = 1/A  on inner boundary
//       du/dn = -1/B^2  on outer boundary
//
//  Exact solution:  u = 1 / sqrt(x^2 + y^2)
//

#define DEBUG_ORDER

#include "hermes2d.h"
#include "solver_umfpack.h"

const double A = 0.05;
const double B = 0.5;


inline double int_w2_u_v(RealFunction* w, RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 2*w->get_fn_order() + fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  w->set_quad_order(o, FN_VAL);
  fu->set_quad_order(o, FN_VAL);
  fv->set_quad_order(o, FN_VAL);

  double *wval = w->get_fn_values();  
  double *vval = fv->get_fn_values();
  double *uval = fu->get_fn_values();

  h1_integrate_expression(wval[i] * wval[i] * uval[i] * vval[i]);
  return result;
}

inline double int_w3_v(RealFunction* w, RealFunction* fu, RefMap* ru)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 3*w->get_fn_order() + fu->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  w->set_quad_order(o, FN_VAL);
  fu->set_quad_order(o, FN_VAL);

  double *wval = w->get_fn_values();  
  double *vval = fu->get_fn_values();

  h1_integrate_expression(wval[i] * wval[i] * wval[i] * vval[i]);
  return result;
}


Solution uprev;

scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv) + int_w2_u_v(&uprev, fu, fv, ru, rv);
}

scalar linear_form_surf(RealFunction* fv, RefMap* rv, EdgePos* ep)
{
  return -1.0/(sqr(B)) * surf_int_v(fv, rv, ep);
}

// scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
// {
//   return int_grad_u_grad_v(fu, fv, ru, rv);
// }
// 
// scalar linear_form(RealFunction* fv, RefMap* rv)
// {
//   return -int_w3_v(&uprev, fv, rv);
// }
//  
// scalar linear_form_surf(RealFunction* fv, RefMap* rv, EdgePos* ep)
// {
//   return -1.0/(sqr(B)) * surf_int_v(fv, rv, ep);
// }



int bc_types(int marker)
{
  if (marker == 1)
    return BC_ESSENTIAL;
  return BC_NATURAL;
}

scalar bc_values(int marker, double x, double y)
{
  if (marker == 1) return 1.0 / A; 
  else return 1.0 / B;
}

scalar exact(double x, double y, scalar& dx, scalar& dy)
{
  return 1.0 / sqrt(sqr(x) + sqr(y));
}


int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("annulus.mesh");
  mesh.refine_all_elements(1);
  mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, 3);
  mesh.refine_all_elements();
//  mesh.refine_all_elements();

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  H1Space space(&mesh, &shapeset);
  //space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(3);
  space.assign_dofs();
  
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form, UNSYM, ANY, 1, &uprev);
  //wf.add_liform(0, linear_form, ANY, 1, &uprev);
  //wf.add_liform_surf(0, linear_form_surf, 2);
  
  UmfpackSolver umfpack;
  LinSystem ls(&wf, &umfpack);
  ls.set_spaces(1, &space);
  ls.set_pss(1, &pss);
  
  ScalarView view("Iteration", 0, 0, 880, 800);
  ScalarView errview("Error", 900, 0, 880, 800);
  
  //uprev.set_exact(&mesh, exact);
  uprev.set_zero(&mesh);

  ExactSolution exsln(&mesh, exact);
  view.show(&exsln);
  view.wait_for_keypress();

  GnuplotGraph graph;
  graph.set_captions("L2 error","iterations","error [%]");
  graph.add_row("fixed point");

  double residuum, error;  
  int it = 0;
  do
  {
    info("------------------ it = %d------------------------------", it++);

    Solution sln;
    ls.assemble();

    ls.solve(1, &sln);
    info("residuum: %g", residuum = l2_error(&uprev, &sln)); 
    info("error: %g", error = l2_error(&sln, &exsln));
    info("");
    
//     graph.add_values(0, it, error);
//     graph.save("conv.txt");

    view.show(&sln);
    DiffFilter err(&sln, &exsln);
    errview.show(&err);
    //errview.wait_for_keypress();

    uprev = sln;
  }
  while (error > 1e-5);
  
  View::wait();
  return 0;
}
