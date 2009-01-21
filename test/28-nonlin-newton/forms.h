#ifndef __FORMS_H
#define __FORMS_H


inline double int_grad_U_grad_v(Solution* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double *dudx, *dudy, *dvdx, *dvdy;
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);

  h1_integrate_dd_expression(dudx[i] * t_dvdx + dudy[i] * t_dvdy);
  return result;
}


inline double int_u3_v(RealFunction* fu, RealFunction* vi, RefMap* ru, RefMap* ri)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 3*fu->get_fn_order() + vi->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  fu->set_quad_order(o, FN_VAL);
  vi->set_quad_order(o, FN_VAL);

  double* uval = fu->get_fn_values();
  double* vval = vi->get_fn_values();
  
  h1_integrate_expression(uval[i]*uval[i]*uval[i]*vval[i]);
  return result;
}


inline double int_u2_vj_vi(RealFunction* fu, RealFunction* vj, RealFunction* vi,
                           RefMap* ru, RefMap* rj, RefMap* ri)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 2*fu->get_fn_order() + vj->get_fn_order() + vi->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  fu->set_quad_order(o);
  vj->set_quad_order(o);
  vi->set_quad_order(o);
  
  double* uval = fu->get_fn_values();
  double* jval = vj->get_fn_values();
  double* ival = vi->get_fn_values();

  h1_integrate_expression(sqr(uval[i]) * jval[i] * ival[i]);
  return result;
}


#endif
