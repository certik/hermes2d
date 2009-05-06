// Neumann BC test: 
//
// Exact solution: u = x^2
//
// PDE: -\Delta u = -2,  \Omega = (0,1)^2
// 
// u = 0      on  x = 0
// du/dn = 0  on  y = 0, y = 1
// du/dn = 2  on  x = 1
//
// Weak formulation: \int_\Omega \nabla u \cdot \nabla v = -2 \int\Omega v + 2 \int_2 v
//

#include "hermes2d.h"
#include "solver_umfpack.h"


int bc_types(int marker)
{
  if (marker == 4) return BC_ESSENTIAL;
  return BC_NATURAL;
}

scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return -2.0 * int_v(fv, rv);
}

scalar linear_form_surf_marker_2(RealFunction* fv, RefMap* rv, EdgePos* ep)
{
  return 2.0 * surf_int_v(fv, rv, ep);
}


int main(int argc, char* argv[])
{
  Mesh mesh;
  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);
  H1Space space(&mesh, &shapeset);
  Solution sln;
  UmfpackSolver umfpack;
  
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form);
  wf.add_liform(0, linear_form);
  wf.add_liform_surf(0, linear_form_surf_marker_2, 2);
  
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);
  
  mesh.load("neumann1.mesh");
  mesh.refine_element(0);
  mesh.refine_all_elements();

  space.set_bc_types(bc_types);
  space.set_uniform_order(3);
  space.assign_dofs();

  sys.assemble();
  sys.solve(1, &sln);
  
  ScalarView view1("Solution - triangles", 200, 150, 1000, 900);
  view1.show(&sln);
  
  mesh.load("neumann2.mesh");
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  space.set_bc_types(bc_types);
  space.set_uniform_order(5);
  space.assign_dofs();

  sys.assemble();
  sys.solve(1, &sln);

  ScalarView view2("Solution - quads", 200, 150, 1000, 900);
  view2.show(&sln);

  View::wait();
  return 0;
}
