// integrals needed in stabilization terms
// maybe later included into integrals_h1.h

inline double int_stab_0_0(RealFunction* w1, RealFunction* w2,
                           RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = fu->get_fn_order() + fv->get_fn_order() + 2 * w1->get_fn_order() + ru->get_inv_ref_order(); 
  limit_order(o);

  w1->set_quad_order(o, FN_VAL);
  w2->set_quad_order(o, FN_VAL);
  fu->set_quad_order(o, FN_ALL);
  fv->set_quad_order(o);

  double *dvdx, *dvdy;
  fv->get_dx_dy_values(dvdx, dvdy);
  double* w1val = w1->get_fn_values();
  double* w2val = w2->get_fn_values();

  double *dudxx, *dudxy, *dudyy;
  dudxx = fu->get_dxx_values();
  dudxy = fu->get_dxy_values();    
  dudyy = fu->get_dyy_values();  
  double *dudx, *dudy;
  fu->get_dx_dy_values(dudx, dudy);
  double* uval = fu->get_fn_values();

  double result = 0.0;
  double3* pt = quad->get_points(o);
  int np = quad->get_num_points(o);
  double2x2 *mu, *mv;
  if (ru->is_jacobian_const()) 
  {
    for (int i = 0; i < np; i++) 
    {
      mu = ru->get_const_inv_ref_map(); 
      mv = rv->get_const_inv_ref_map(); 
      double a = (sqr((*mu)[0][0]) + sqr((*mu)[1][0]));
      double b = (sqr((*mu)[0][1]) + sqr((*mu)[1][1]));
      double c = 2.0 * ((*mu)[0][0]*(*mu)[0][1] + (*mu)[1][0]*(*mu)[1][1]);
      result += pt[i][2] * ( ( - ((dudxx[i]*a + dudxy[i]*c + dudyy[i]*b ) / Re)
                               + (w1val[i] * t_dudx + w2val[i] * t_dudy)
                               + (uval[i] / tau) 
                             )
                             * ((w1val[i] * t_dvdx + w2val[i] * t_dvdy)) ); 
    }
    result *= ru->get_const_jacobian();
  }
  else
  {
    mu = ru->get_inv_ref_map(o); 
    mv = rv->get_inv_ref_map(o); 
    double3x2 *mm;
    mm = ru->get_second_ref_map(o); 
    double* jac = ru->get_jacobian(o); 
    for (int i = 0; i < np; i++, mu++, mv++, mm++) 
    {
      double a = (sqr((*mu)[0][0]) + sqr((*mu)[1][0]));
      double b = (sqr((*mu)[0][1]) + sqr((*mu)[1][1]));
      double c = 2.0 * ((*mu)[0][0]*(*mu)[0][1] + (*mu)[1][0]*(*mu)[1][1]);
      double coefx = (*mm)[0][0] + (*mm)[2][0];
      double coefy = (*mm)[0][1] + (*mm)[2][1];
      result += pt[i][2] * (jac[i]) * ( ( - (( dudx[i]*coefx + dudy[i]*coefy + dudxx[i]*a + dudxy[i]*c + dudyy[i]*b ) / Re)
                                          + (w1val[i] * t_dudx + w2val[i] * t_dudy)
                                          + (uval[i] / tau) 
                                        )
                                        * ((w1val[i] * t_dvdx + w2val[i] * t_dvdy)) ); 
    }
  }
  return result;
}


inline double int_dudx_w_nabla_v(RealFunction* w1, RealFunction* w2,
                              RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = fu->get_fn_order() + fv->get_fn_order() + w1->get_fn_order() + ru->get_inv_ref_order(); 
  limit_order(o);

  w1->set_quad_order(o, FN_VAL);
  w2->set_quad_order(o, FN_VAL);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double *dudx, *dudy;
  fu->get_dx_dy_values(dudx, dudy);
  double *dvdx, *dvdy;
  fv->get_dx_dy_values(dvdx, dvdy);
  double* w1val = w1->get_fn_values();
  double* w2val = w2->get_fn_values();

  h1_integrate_dd_expression((t_dudx) * (w1val[i] * t_dvdx + w2val[i] * t_dvdy));
  return result;
}

inline double int_dudy_w_nabla_v(RealFunction* w1, RealFunction* w2,
                              RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  
  int o = fu->get_fn_order() + fv->get_fn_order() + w1->get_fn_order() + ru->get_inv_ref_order(); 
  limit_order(o);

  w1->set_quad_order(o, FN_VAL);
  w2->set_quad_order(o, FN_VAL);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double *dudx, *dudy;
  fu->get_dx_dy_values(dudx, dudy);
  double *dvdx, *dvdy;
  fv->get_dx_dy_values(dvdx, dvdy);
  double* w1val = w1->get_fn_values();
  double* w2val = w2->get_fn_values();

  h1_integrate_dd_expression((t_dudy) * (w1val[i] * t_dvdx + w2val[i] * t_dvdy));
  return result;
}


