// Newton BC test:
//
// Exact solution: u = x^2
//
// PDE: -\Delta u = -2,  \Omega = (0,1)^2
//
// u = 0          on  x = 0
// du/dn = 0      on  y = 0, y = 1
// u + du/dn = 3  on  x = 1
//
// Weak formulation:  \int_\Omega \nabla u \cdot \nabla v + \int_2 uv = -2\int\Omega v + 3\int_2 v
//

#include "hermes2d.h"


int bc_types(int marker)
{
  return (marker == 4) ? BC_ESSENTIAL : BC_NATURAL;
}

scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar bilinear_form_surf(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv, EdgePos* ep)
{
  return (ep->marker == 2) ? surf_int_u_v(fu, fv, ru, rv, ep) : 0.0;
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return -2.0 * int_v(fv, rv);
}

scalar linear_form_surf(RealFunction* fv, RefMap* rv, EdgePos* ep)
{
  return (ep->marker == 2) ? 3.0 * surf_int_v(fv, rv, ep) : 0.0;
}


int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);

  Mesh mesh;
  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);
  H1Space space(&mesh, &shapeset);
  Solution sln;
  
  DiscreteProblem dp;
  dp.set_num_equations(1);
  dp.set_spaces(1, &space);
  dp.set_pss(1, &pss);
  dp.set_bilinear_form(0, 0, NULL, bilinear_form, bilinear_form_surf);
  dp.set_linear_form(0, linear_form, linear_form_surf);
  
  mesh.load("newton1.mesh");
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  space.set_bc_types(bc_types);
  space.set_uniform_order(5);
  space.assign_dofs();

  dp.create_matrix();
  dp.assemble_matrix_and_rhs();
  dp.solve_system(1, &sln);

  ScalarView view1("Solution - triangles", 200, 150, 1000, 900);
  view1.show(&sln);
  
  mesh.load("newton2.mesh");
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  space.set_bc_types(bc_types);
  space.set_uniform_order(5);
  space.assign_dofs();

  dp.create_matrix();
  dp.assemble_matrix_and_rhs();
  dp.solve_system(1, &sln);

  ScalarView view2("Solution - quads", 200, 150, 1000, 900);
  view2.show(&sln);

  hermes2d_finalize();
  return 0;
}
