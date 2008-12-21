#include "hermes2d.h"

const double Re = 300;
const double tau = 0.05;


int xvel_bc_type(int marker)
  { return (marker != 11) ? BC_ESSENTIAL : BC_NONE; }

scalar xvel_bc_value(int marker, double x, double y)
  { return (marker > 3) ? 1 : 0; }

int yvel_bc_type(int marker)
  { return (marker != 11) ? BC_ESSENTIAL : BC_NONE; }

int press_bc_type(int marker)
  { return BC_NONE; }


Solution xprev;
Solution yprev;


scalar bilinear_form_both_0_0_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv) / Re + 
           int_u_v(fu, fv, ru, rv) / tau + 
           int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv); }

scalar bilinear_form_sym_0_0_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv) / Re + 
           int_u_v(fu, fv, ru, rv) / tau; }

scalar bilinear_form_unsym_0_0_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv); }

scalar bilinear_form_unsym_0_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return -int_u_dvdx(fu, fv, ru, rv); }

scalar bilinear_form_unsym_1_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return -int_u_dvdy(fu, fv, ru, rv); }

scalar bilinear_form_unsym_2_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_dudx_v(fu, fv, ru, rv); }

scalar bilinear_form_unsym_2_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_dudy_v(fu, fv, ru, rv); }

scalar linear_form_0(RealFunction* fv, RefMap* rv)
  { return int_u_v(&xprev, fv, rv, rv) / tau; }

scalar linear_form_1(RealFunction* fv, RefMap* rv)
  { return int_u_v(&yprev, fv, rv, rv) / tau; }


int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);
  
  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);
  
  Mesh mesh;
  mesh.load("airfoil.mesh");
  //mesh.refine_all_elements();
  
  int ndofs = 0;
  
  H1Space xvel(&mesh, &shapeset);
  xvel.set_bc_types(xvel_bc_type);
  xvel.set_bc_values(xvel_bc_value);
  xvel.set_uniform_order(2);
  ndofs += xvel.assign_dofs(ndofs);
  
  H1Space yvel(&mesh, &shapeset);
  yvel.set_bc_types(yvel_bc_type);
  yvel.set_uniform_order(2);
  ndofs += yvel.assign_dofs(ndofs);
  
  H1Space press(&mesh, &shapeset);
  press.set_bc_types(press_bc_type);
  press.set_uniform_order(1);
  ndofs += press.assign_dofs(ndofs);

  xprev.set_zero(&mesh);
  yprev.set_zero(&mesh);
  
  DiscreteProblem dp;
  dp.set_num_equations(3);
  dp.set_spaces(3, &xvel, &yvel, &press);
  dp.set_pss(1, &pss);
  dp.set_external_fns(2, &xprev, &yprev);

  dp.set_bilinear_form(0, 0, bilinear_form_unsym_0_0_1_1, bilinear_form_sym_0_0_1_1);
  dp.set_bilinear_form(1, 1, bilinear_form_unsym_0_0_1_1, bilinear_form_sym_0_0_1_1);
  //dp.set_bilinear_form(0, 0, bilinear_form_both_0_0_1_1);
  //dp.set_bilinear_form(1, 1, bilinear_form_both_0_0_1_1);
  
  dp.set_bilinear_form(1, 2, bilinear_form_unsym_1_2);
  dp.set_bilinear_form(0, 2, bilinear_form_unsym_0_2);
  
  //dp.set_bilinear_form(2, 0, bilinear_form_unsym_2_0);
  //dp.set_bilinear_form(2, 1, bilinear_form_unsym_2_1);
  dp.set_bilinear_form(2, 0, BF_ANTISYM);
  dp.set_bilinear_form(2, 1, BF_ANTISYM);
  
  dp.set_linear_form(0, linear_form_0);
  dp.set_linear_form(1, linear_form_1);
  
  dp.create_matrix();
  
  //MatrixView mv;
  VectorView vview("Velocity");
  ScalarView pview("Pressure");
  
  for (int i = 0; i < 10; i++)
  {
    Solution psln;
    dp.assemble_matrix_and_rhs();
    dp.solve_system(3, &xprev, &yprev, &psln);
        
    //mv.show(&ep);
    //mv.wait_for_close();
    
    vview.show(&xprev, &yprev);
    pview.show(&psln);
    
    printf("\n");
  }
  
  hermes2d_finalize();
  return 0;
}


// gcc asm 0.81s, total (10 iterations) 21.4s

/*scalar test_vectorize(int np, double* weights[3], double* dudx, double* dudy, double* dvdx, double* dvdy, double* m[4])
{
  scalar result = 0.0;
  for (int i = 0; i < np; i++)
  {
    result += weights[2][i] * ((dudx[i]*m[0][i] + dudy[i]*m[1][i]) * (dvdx[i]*m[0][i] + dudy[i]*m[1][i]) +
                               (dudx[i]*m[2][i] + dudy[i]*m[3][i]) * (dvdx[i]*m[2][i] + dudy[i]*m[3][i]));
    //result += weights[i] * dudx[i]*m1[i][0];
  }
  return result;
}*/
