/********** Definition of Jacobian matrices and residual vectors ***********/ // Residuum for the Euler time discretization

// Residuum for the Euler time discretization
inline double F_euler(RealFunction* Tprev, RealFunction* Titer, RealFunction* fu, RefMap* ru)
{
  Quad2D* quad = fu->get_quad_2d();
  RefMap* rv = ru;

  int o = 3 * Titer->get_fn_order() + fu->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Tprev->set_quad_order(o, FN_VAL);
  Titer->set_quad_order(o);
  fu->set_quad_order(o);

  double* Titer_val = Titer->get_fn_values();
  double* Tprev_val = Tprev->get_fn_values();
  double* uval = fu->get_fn_values();
  double *dTiter_dx, *dTiter_dy, *dudx, *dudy;
  Titer->get_dx_dy_values(dTiter_dx, dTiter_dy);
  fu->get_dx_dy_values(dudx, dudy);

  // u is a test function
  double result = 0.0;
  h1_integrate_dd_expression(( HEATCAP*(Titer_val[i] - Tprev_val[i])*uval[i]/TAU +
                               lam(Titer_val[i]) * (dTiter_dx[i]*t_dudx + dTiter_dy[i]*t_dudy)));

  return result;
}

// Jacobian matrix for the implicit Euler time discretization
inline double J_euler(RealFunction* Titer, RealFunction* fu,
                RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 2 * Titer->get_fn_order() + fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Titer->set_quad_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double* Titer_val = Titer->get_fn_values();
  double* uval = fu->get_fn_values();
  double* vval = fv->get_fn_values();

  double *dTiter_dx, *dTiter_dy, *dudx, *dudy, *dvdx, *dvdy;
  Titer->get_dx_dy_values(dTiter_dx, dTiter_dy);
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);

  // u is a basis function, v a test function
  double result = 0.0;
  h1_integrate_dd_expression(( HEATCAP * uval[i] * vval[i] / TAU +
                               dlam_dT(Titer_val[i]) * uval[i] * (dTiter_dx[i]*t_dvdx + dTiter_dy[i]*t_dvdy) +
                               lam(Titer_val[i]) * (t_dudx*t_dvdx + t_dudy*t_dvdy)));

  return result;
}

// Residuum for the Crank-Nicolson time discretization
inline double F_cranic(RealFunction* Tprev, RealFunction* Titer, RealFunction* fu, RefMap* ru)
{
  Quad2D* quad = fu->get_quad_2d();
  RefMap* rv = ru;

  int o = 3 * Titer->get_fn_order() + fu->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Tprev->set_quad_order(o);
  Titer->set_quad_order(o);
  fu->set_quad_order(o);

  double* Titer_val = Titer->get_fn_values();
  double* Tprev_val = Tprev->get_fn_values();
  double* uval = fu->get_fn_values();
  double *dTiter_dx, *dTiter_dy, *dTprev_dx, *dTprev_dy, *dudx, *dudy;
  Titer->get_dx_dy_values(dTiter_dx, dTiter_dy);
  Tprev->get_dx_dy_values(dTprev_dx, dTprev_dy);
  fu->get_dx_dy_values(dudx, dudy);

  // u is a test function
  double result = 0.0;
  h1_integrate_dd_expression(( HEATCAP * (Titer_val[i] - Tprev_val[i]) * uval[i] / TAU +
                               0.5 * lam(Titer_val[i]) * (dTiter_dx[i]*t_dudx + dTiter_dy[i]*t_dudy) +
                               0.5 * lam(Tprev_val[i]) * (dTprev_dx[i]*t_dudx + dTprev_dy[i]*t_dudy)
                            ));

  return result;
}

// Jacobian matrix for the Crank-Nicolson time discretization
inline double J_cranic(RealFunction* Titer, RealFunction* fu,
                RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 2 * Titer->get_fn_order() + fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Titer->set_quad_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double* Titer_val = Titer->get_fn_values();
  double* uval = fu->get_fn_values();
  double* vval = fv->get_fn_values();

  double *dTiter_dx, *dTiter_dy, *dudx, *dudy, *dvdx, *dvdy;
  Titer->get_dx_dy_values(dTiter_dx, dTiter_dy);
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);

  // u is a basis function, v a test function
  double result = 0.0;
  h1_integrate_dd_expression(( HEATCAP * uval[i] * vval[i] / TAU +
                               0.5 * dlam_dT(Titer_val[i]) * uval[i] * (dTiter_dx[i]*t_dvdx + dTiter_dy[i]*t_dvdy) +
                               0.5 * lam(Titer_val[i]) * (t_dudx*t_dvdx + t_dudy*t_dvdy)
                            ));

  return result;
}
