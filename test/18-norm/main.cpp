#include "hermes2d.h"


double exact(double x, double y, double& dx, double& dy)
{
  dx = 2*x;
  dy = 2*y;
  return sqr(x) + sqr(y);
}


double bc_values(int marker, double x, double y)
{
  return sqr(x) + sqr(y);
}

scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return -4 * int_v(fv, rv);
}


int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);

  Mesh mesh;
  mesh.load("square1.mesh");
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
  space.set_bc_values(bc_values);
  space.set_uniform_order(4);
  space.assign_dofs();
  
  DiscreteProblem dp;
  dp.set_num_equations(1);
  dp.set_spaces(1, &space);
  dp.set_pss(1, &pss);
  dp.set_bilinear_form(0, 0, bilinear_form);
  dp.set_linear_form(0, linear_form);

  Solution sln;
  dp.create_matrix();
  dp.assemble_matrix_and_rhs();
  dp.solve_system(1, &sln);

  ScalarView view("Solution", 200, 150, 1000, 900);
  view.show(&sln);
  
  printf("%.18g\n", h1_error_norm(&sln, exact));

  hermes2d_finalize();
  return 0;
}

