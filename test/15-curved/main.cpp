//
//  Curved element test
//
//  Domain:  (x^2 + y^2 < 1)  &  (x > 0)  &  (y > 0)
//
//  Exact solution:  u = x^2
//
//  PDE: -\Delta u = -2
//
//  BC:          u = 0 on (x = 0)
//           du/dn = 0 on (y = 0)
//           du/dn = 2x^2 / sqrt(x^2 + y^2) on (x^2 + y^2 = 1)
//

#include "hermes2d.h"
#include "solver_umfpack.h"


int bc_types(int marker)
{
  return (marker == 3) ? BC_ESSENTIAL : BC_NATURAL;
}

double bc_values(int marker, double x, double y)
{
  return (marker == 2) ? 2.0*x*x / sqrt(x*x + y*y) : 0.0;
}

scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return -2.0 * int_v(fv, rv);
}

scalar linear_form_surf(RealFunction* fv, RefMap* rv, EdgePos* ep)
{
  return surf_int_G_v(fv, rv, ep);
}


int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("domain.mesh");
  
  Mesh mesh2;
  mesh2.copy(&mesh);
  mesh2.refine_all_elements();
  Solution dummy;
  dummy.set_zero(&mesh2);

  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(6);
  space.assign_dofs();
  
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form, UNSYM, ANY, 1, &dummy);
  wf.add_liform(0, linear_form);
  wf.add_liform_surf(0, linear_form_surf, ANY, 1, &dummy);
  
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);
  
  Solution sln;
  sys.assemble();
  sys.solve(1, &sln);

  ScalarView view1("Solution 1");
  view1.show(&sln);
  
  View::wait();
  return 0;
}
