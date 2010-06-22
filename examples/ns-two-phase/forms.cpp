

template<typename Real, typename Scalar>
double int_sign_u_v(double a, double b, Func<Real>* l, Func<Real>* fu, Func<Real>* fv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 5;
  limit_order(o);

  l->set_quad_order(o, H2D_FN_VAL); 
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double* uval = fu->get_fn_values();
  double* vval = fv->get_fn_values();
  double* lval = l->get_fn_values();

  // FIXME
  //h1_integrate_expression(((lval[i] < 0.0) ? b : a) * uval[i] * vval[i]);
  double result = 0;
  return result;
}

template<typename Real, typename Scalar>
double int_sign_grad_u_grad_v(double a, double b, Func<Real>* l, Func<Real>* fu, Func<Real>* fv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 5;
  limit_order(o);
  l->set_quad_order(o, H2D_FN_VAL);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double *dudx, *dudy, *dvdx, *dvdy;
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);
  double* lval = l->get_fn_values();

  // FIXME: Where is number i in lval[i]?
  //h1_integrate_dd_expression(((lval[i] < 0.0) ? b : a) * (t_dudx * t_dvdx + t_dudy * t_dvdy));
  double result = 0;
  return result;
}

template<typename Real, typename Scalar>
double int_sign_w_nabla_u_v(double a, double b, Func<Real>* l, Func<Real>* w1, Func<Real>* w2,
			    Func<Real>* fu, Func<Real>* fv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 5;
  limit_order(o);

  l->set_quad_order(o, H2D_FN_VAL);
  w1->set_quad_order(o, H2D_FN_VAL);
  w2->set_quad_order(o, H2D_FN_VAL);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double *dudx, *dudy;
  fu->get_dx_dy_values(dudx, dudy);
  double* lval = l->get_fn_values();
  double* vval = fv->get_fn_values();
  double* w1val = w1->get_fn_values();
  double* w2val = w2->get_fn_values();

  // FIXME: the same as above, could not find the define of number 'i'.
  //h1_integrate_dd_expression(((lval[i] < 0.0) ? b : a) * (w1val[i] * t_dudx + w2val[i] * t_dudy) * vval[i]);
  double result = 0;
  return result;
}

template<typename Real, typename Scalar>
double int_sign_v(double a, double b, Func<Real>* l, Func<Real>* fu)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 5;
  limit_order(o);
 
  l->set_quad_order(o, H2D_FN_VAL);  
  fu->set_quad_order(o);

  double* uval = fu->get_fn_values(); 
  double* lval = l->get_fn_values();

  // FIXME: the same as above, could not find the define of number 'i'.
  //h1_integrate_expression(((lval[i] < 0.0) ? b : a) * uval[i]);)
  double result = 0;
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_0_0(int n, double *wt, Func<Real>* fu, Func<Real>* fv, 
                         Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  Func<Scalar>* lprev = ext->fn[0];
  Func<Scalar>* xprev = ext->fn[1];
  Func<Scalar>* yprev = ext->fn[2];

  return int_sign_grad_u_grad_v(Nu1, Nu2, lprev, fu, fv) + int_sign_u_v(Ro1, Ro2, lprev, fu, fv)/tau + 
         int_sign_w_nabla_u_v( Ro1, Ro2, lprev, xprev, yprev, fu, fv); 
}

// FIXME: Please check the following functions.
//scalar bilinear_form_0_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
//  { return int_dudx_v(fu, fv, ru ,rv ); }
template<typename Real, typename Scalar>
Scalar bilinear_form_0_2(int n, double *wt, Func<Real>* fu, Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext)
  { return int_dudx_v<Real, Scalar>(n, wt, fu, fv); }
 
// FIXME: Please check the following functions.
//scalar bilinear_form_1_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
//  { return int_dudy_v(fu, fv, ru ,rv ); } 
template<typename Real, typename Scalar>
Scalar bilinear_form_1_2(int n, double *wt, Func<Real>* fu, Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext)
  { return int_dudy_v<Real, Scalar>(n, wt, fu, fv); }

// FIXME: 
//scalar bilinear_form_2_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv) 
//  { return int_dudx_v(fu, fv, ru, rv);}
template<typename Real, typename Scalar>
Scalar bilinear_form_2_0(int n, double *wt, Func<Real>* fu, Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext) 
  { return int_dudx_v<Real, Scalar>(n, wt, fu, fv);}

// FIXME:
//scalar bilinear_form_2_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)  
//  { return int_dudy_v(fu, fv, ru, rv); }
template<typename Real, typename Scalar>
Scalar bilinear_form_2_1(int n, double *wt, Func<Real>* fu, Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext)  
  { return int_dudy_v<Real, Scalar>(n, wt, fu, fv); }

// FIXME:
//scalar bilinear_form_3_3(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
//  { return  int_u_v(fu, fv, ru, rv)/tau + int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv); }
template<typename Real, typename Scalar>
Scalar bilinear_form_3_3(int n, double *wt, Func<Real>* fu, Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  Func<Scalar>* xprev = ext->fn[0];
  Func<Scalar>* yprev = ext->fn[1];

  return  int_u_v<Real, Scalar>(n, wt, fu, fv)/tau + int_w_nabla_u_v(n, wt, &xprev, &yprev, fu, fv); }

// FIXME:
//scalar linear_form_0(RealFunction* fv, RefMap* rv)
//  { return int_sign_u_v(Ro1, Ro2, &lprev, &xprev, fv, rv, rv)/tau; }
template<typename Real, typename Scalar>
Scalar linear_form_0(int n, double *wt, Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  Func<Scalar>* lprev = ext->fn[0];
  Func<Scalar>* xprev = ext->fn[1];

  return int_sign_u_v(Ro1, Ro2, lprev, xprev, fv)/tau; }

// FIXME:
//scalar linear_form_1(RealFunction* fv, RefMap* rv)
//  { return int_sign_u_v(Ro1, Ro2, &lprev, &yprev, fv, rv, rv)/tau - int_sign_v(10*Ro1, 10*Ro2, &lprev, fv, rv); }
template<typename Real, typename Scalar>
Scalar linear_form_1(int n, double *wt, Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  Func<Scalar>* lprev = ext->fn[0];
  Func<Scalar>* yprev = ext->fn[1];

  return int_sign_u_v(Ro1, Ro2, lprev, &yprev, fv)/tau - int_sign_v(10*Ro1, 10*Ro2, lprev, fv); 
}

// FIXME:
//Scalar linear_form_3(RealFunction* fv, RefMap* rv)
//  { return int_u_v(&lprev, fv, rv, rv)/tau; }
template<typename Real, typename Scalar>
Scalar linear_form_3(int n, double *wt, Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext)
{  
  Func<Scalar>* lprev = ext->fn[0];

  return int_u_v<Real, Scalar>(n, wt, lprev, fv)/tau; 
}
