#include "hermes2d.h"

const double tau = 0.05;


int temp_bc_types(int marker)
  { return (marker != 2) ? BC_ESSENTIAL : BC_NATURAL; }

scalar temp_bc_values(int marker, double x, double y)
  { return (marker == 1 || marker == 3) ? 1.0 : 0.0; }


Solution tprev;

scalar bilinear_form_unsym(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_u_v(fu, fv, ru, rv) / tau +
         0.02 * int_grad_u_grad_v(fu, fv, ru, rv) +
         int_dudx_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return int_u_v(&tprev, fv, rv, rv) / tau;
}



int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);
  
  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);
  
  Mesh mesh;
  mesh.load("channel.mesh");
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  //mesh.refine_all_elements();
  //mesh.refine_all_elements();
  
  H1Space temp(&mesh, &shapeset);
  temp.set_bc_types(temp_bc_types);
  temp.set_bc_values(temp_bc_values);
  temp.set_uniform_order(3);
  int ndofs = temp.assign_dofs();
  
  tprev.set_zero(&mesh);
  
  DiscreteProblem dp;
  dp.set_num_equations(1);
  dp.set_spaces(1, &temp);
  dp.set_pss(1, &pss);
  dp.set_external_fns(1, &tprev);
  dp.set_bilinear_form(0, 0, bilinear_form_unsym);
  dp.set_linear_form(0, linear_form);
  dp.create_matrix();
  
  ScalarView tview("Temperature", 100, 300, 1400, 600);
  for (int i = 0; i < 100; i++)
  {
    dp.assemble_matrix_and_rhs();
    dp.solve_system(1, &tprev);
    
    tview.show(&tprev);
    printf("\n\n");
  }
  
  hermes2d_finalize();
  return 0;
}
