#include "hermes2d.h"


scalar bilinear_form_unsym(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv); }

scalar linear_form(RealFunction* fv, RefMap* rv)
  { 
    return int_v(fv, rv); 
  }


int main(int argc, char* argv[])
{
  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);
  
  Mesh mesh0, mesh1;
  mesh0.load("square1.mesh");
  mesh1.copy(&mesh0);
  mesh1.refine_towards_vertex(3, 10);
  
  Solution dummy;
  dummy.set_zero(&mesh1);

  H1Space space(&mesh0, &shapeset);
  space.set_uniform_order(2);
  space.assign_dofs();

  DiscreteProblem dp;
  dp.set_num_equations(1);
  dp.set_spaces(1, &space);
  dp.set_pss(1, &pss);
  dp.set_external_fns(1, &dummy);

  dp.set_bilinear_form(0, 0, bilinear_form_unsym);
  dp.set_linear_form(0, linear_form);
  
  Solution sln;
  dp.create_matrix();
  dp.assemble_matrix_and_rhs();
  dp.save_rhs_matlab("bad.m", "rhs");
  dp.solve_system(1, &sln);

  ScalarView view;
  view.show(&sln);

  View::wait();
  return 0;
}
