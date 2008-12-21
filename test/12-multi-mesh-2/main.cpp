#include "hermes2d.h"
#include <unistd.h>


const double Re = 500;
const double tau = 0.05;


int xvel_bc_types(int marker)
  { return (marker != 4) ? BC_ESSENTIAL : BC_NONE; }
  
scalar xvel_bc_values(int marker, double x, double y)
  { return (marker == 8) ? 1.5*(1.0 - pow(2*y-1, 20)) : 0; }

int yvel_bc_types(int marker)
  { return (marker != 4) ? BC_ESSENTIAL : BC_NONE; }

int press_bc_types(int marker)
  { return BC_NONE; }

int temp_bc_types(int marker)
  { return (marker == 8) ? BC_ESSENTIAL : BC_NATURAL; }
  //{ return (marker == 8 || marker == 2 || marker == 6) ? BC_ESSENTIAL : BC_NATURAL; }

//scalar temp_bc_values(int marker, double x, double y)
//  { return (marker == 2 || marker == 6) ? 1.0 : 0.0; }
  

  

Solution xprev;
Solution yprev;
Solution tprev;

scalar bilinear_form_unsym_0_0_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv) / Re + 
           int_u_v(fu, fv, ru, rv) / tau + 
           int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv); }

scalar bilinear_form_unsym_0_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return -int_u_dvdx(fu, fv, ru, rv); }

scalar bilinear_form_unsym_1_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return -int_u_dvdy(fu, fv, ru, rv); }

scalar bilinear_form_unsym_2_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_dudx_v(fu, fv, ru, rv); }

scalar bilinear_form_unsym_2_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_dudy_v(fu, fv, ru, rv); }
  
scalar bilinear_form_unsym_3_3(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return 0.03 * int_grad_u_grad_v(fu, fv, ru, rv) +
           int_u_v(fu, fv, ru, rv) / tau +
           int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv); }
  
scalar linear_form_0(RealFunction* fv, RefMap* rv)
  { return int_u_v(&xprev, fv, rv, rv /*fixme*/) / tau; }

scalar linear_form_1(RealFunction* fv, RefMap* rv)
  { return int_u_v(&yprev, fv, rv, rv /*fixme*/) / tau; }
  
scalar linear_form_3(RealFunction* fv, RefMap* rv)
  { return int_u_v(&tprev, fv, rv, rv /*fixme*/) / tau; }
  
scalar linear_form_3_surf(RealFunction* fv, RefMap* refmap, EdgePos* ep)
  { return (ep->marker == 2 || ep->marker == 6) ? 2*surf_int_v(fv, refmap, ep) : 0; }



int ref_crit(Element* e)
{ 
  if (e->id < 5) return 0;
  return 1;
}

int ref_crit_x1(Element* e)
  { return (e->en[0]->marker > 0 || e->en[2]->marker > 0) ? 1 : -1; }

int ref_crit_y1(Element* e)
  { return (e->en[0]->marker > 0 || e->en[2]->marker > 0) ? 1 : -1; }

int ref_crit_p1(Element* e)
  { return (e->en[0]->marker == 2 || e->en[2]->marker == 6) ? 1 : -1; }

int ref_crit_p2(Element* e)
  { return ((e->en[0]->marker == 1 || e->en[2]->marker == 7) && (e->vn[0]->x <= 1)) ? 1 : -1; }

int ref_crit_p3(Element* e)
  { return ((e->en[0]->marker == 1 || e->en[2]->marker == 7) && (e->vn[0]->x > 1)) ? 1 : -1; }


int ref_crit_t0(Element* e)
  { return (e->id > 0 && e->id < 6) ? 0 : -1; }

int ref_crit_t1(Element* e)
  { return (e->en[0]->marker == 2 || e->en[2]->marker == 6) ? 1 : -1; }

int ref_crit_t2(Element* e)
  { return (e->en[0]->marker == 3 || e->en[2]->marker == 5) && (e->vn[0]->x < 5) ? 1 : -1; }

int ref_crit_t3(Element* e)
  { return (e->en[0]->marker == 1 || e->en[2]->marker == 7) && (e->vn[0]->x > 1) ? 1 : -1; }

    
int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);
  
  H1ShapesetOrtho shapeset;
  PrecalcShapeset xpss(&shapeset);
  PrecalcShapeset ypss(&shapeset);
  PrecalcShapeset ppss(&shapeset);
  PrecalcShapeset tpss(&shapeset);
  
  /////////////////////////////////////////////////////////////////////////////////////////////
  
  const char* basename = "channel2.mesh";
  
  Mesh xmesh;  
  xmesh.load(basename);
  /*xmesh.refine_element(1, 2);
  xmesh.refine_element(2, 2);
  xmesh.refine_element(3, 2);
  xmesh.refine_element(4, 2);*/
  //xmesh.refine_all_elements(2);
  //xmesh.refine_all_elements();
  xmesh.refine_by_criterion(ref_crit, 1);
  xmesh.refine_by_criterion(ref_crit_x1, 8);

  Mesh ymesh;
  ymesh.load(basename);
  ymesh.refine_all_elements();
  //ymesh.refine_all_elements(2);
  //ymesh.refine_by_criterion(ref_crit, 1);
  ymesh.refine_by_criterion(ref_crit_y1, 4);
  
  Mesh pmesh;
  pmesh.load(basename);
  //pmesh.refine_all_elements();
  pmesh.refine_by_criterion(ref_crit, 1);
  pmesh.refine_all_elements(1);
  //pmesh.refine_all_elements(1);  
  pmesh.refine_by_criterion(ref_crit_p1, 2);
  pmesh.refine_by_criterion(ref_crit_p2, 2);
  pmesh.refine_by_criterion(ref_crit_p3, 3);
  
  
  Mesh tmesh;
  tmesh.load(basename);
  //tmesh.refine_all_elements();
  //tmesh.refine_all_elements();
  //tmesh.refine_all_elements();
  tmesh.refine_by_criterion(ref_crit_t0, 1);
  tmesh.refine_by_criterion(ref_crit_t1, 5); 
  tmesh.refine_by_criterion(ref_crit_t2, 3); 
  tmesh.refine_by_criterion(ref_crit_t3, 3); 
  
  /*MeshView mv("Pressure",    0, 586, 1600, 267);
  mv.show(&pmesh);*/
  
  /*xmesh.refine_all_elements();
  ymesh.refine_all_elements();
  pmesh.refine_all_elements();
  tmesh.refine_all_elements();*/
    
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  
  int ndofs = 0;
  
  H1Space xvel(&xmesh, &shapeset);
  xvel.set_bc_types(xvel_bc_types);
  xvel.set_bc_values(xvel_bc_values);
  xvel.set_uniform_order(3);
  ndofs += xvel.assign_dofs(ndofs);
  
  H1Space yvel(&ymesh, &shapeset);
  yvel.set_bc_types(yvel_bc_types);
  yvel.set_uniform_order(3);
  ndofs += yvel.assign_dofs(ndofs);
  
  H1Space press(&pmesh, &shapeset);
  press.set_bc_types(press_bc_types);
  press.set_uniform_order(2);
  ndofs += press.assign_dofs(ndofs);
  
  H1Space temp(&tmesh, &shapeset);
  temp.set_bc_types(temp_bc_types);
  //temp.set_bc_values(temp_bc_values);
  temp.set_uniform_order(3);
  ndofs += temp.assign_dofs(ndofs);
  
  printf("%d\n", ndofs);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  xprev.set_zero(&xmesh);
  yprev.set_zero(&ymesh);
  tprev.set_zero(&tmesh);
  
  DiscreteProblem dp;
  dp.set_num_equations(4);
  dp.set_spaces(4, &xvel, &yvel, &press, &temp);
  dp.set_pss(4, &xpss, &ypss, &ppss, &tpss);
  dp.set_external_fns(3, &xprev, &yprev, &tprev);

  dp.set_bilinear_form(0, 0, bilinear_form_unsym_0_0_1_1);
  dp.set_bilinear_form(0, 2, bilinear_form_unsym_0_2);
  dp.set_bilinear_form(1, 1, bilinear_form_unsym_0_0_1_1);
  dp.set_bilinear_form(1, 2, bilinear_form_unsym_1_2);
  dp.set_bilinear_form(2, 0, bilinear_form_unsym_2_0);
  dp.set_bilinear_form(2, 1, bilinear_form_unsym_2_1);
  dp.set_bilinear_form(3, 3, bilinear_form_unsym_3_3);
  dp.set_linear_form(0, linear_form_0);
  dp.set_linear_form(1, linear_form_1);
  dp.set_linear_form(3, linear_form_3, linear_form_3_surf);
  
  dp.create_matrix();

  ScalarView xview("X velocity",  0, 0,   1600, 267);
  ScalarView yview("Y velocity",  0, 293, 1600, 267);
  ScalarView pview("Pressure",    0, 586, 1600, 267);
  ScalarView tview("Temperature", 0, 879, 1600, 267);
  VectorView vview("Velocity", 0, 0, 1600, 500);

  //xview.show(&xprev, 3);
  //yview.show(&yprev, 3);
  //pview.show(&psln,  3);
  //tview.show(&tprev, 3);

  Solution psln;
  for (int i = 0; i < 100; i++)
  {
    printf("i=%d\n", i);
    dp.assemble_matrix_and_rhs();
    dp.solve_system(4, &xprev, &yprev, &psln, &tprev);

//    if (!(i % 5))
    {
      xview.show(&xprev);
      yview.show(&yprev);
      //vview.show(&xprev, &yprev);
      pview.show(&psln);
      tview.show(&tprev);
    }

    printf("\n\n");
  }
  
  //xview.show(&xprev, 3);
  //yview.show(&yprev, 3);
  //vview.show(&xprev, FN_VAL_0, &yprev, FN_VAL_0);
  //pview.show(&psln, 3);
  //tview.show(&tprev, 3);
  
  hermes2d_finalize();
  return 0;
}
