// Residuum for the implicit Euler time discretization
inline complex F_euler(ScalarFunction* Psi_prev, ScalarFunction* Psi_iter, RealFunction* fu, RefMap* ru)
{
  scalar ii = complex(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Quad2D* quad = fu->get_quad_2d();
  RefMap* rv = ru;

  //not sure here:
  int o = 3 * Psi_iter->get_fn_order() + fu->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Psi_prev->set_quad_order(o, FN_VAL);
  Psi_iter->set_quad_order(o);
  fu->set_quad_order(o);

  scalar* Psi_iter_val = Psi_iter->get_fn_values();
  scalar* Psi_prev_val = Psi_prev->get_fn_values();
  scalar* dPsi_iter_dx, *dPsi_iter_dy;
  double* uval = fu->get_fn_values();
  double *dudx, *dudy;
  Psi_iter->get_dx_dy_values(dPsi_iter_dx, dPsi_iter_dy);
  fu->get_dx_dy_values(dudx, dudy);

  // obtain physical coordinates of int. points
  double* x = ru->get_phys_x(o);
  double* y = ru->get_phys_y(o);

  // u is a test function
  scalar result = 0.0;
  h1_integrate_dd_expression((
    ii * H * (Psi_iter_val[i] - Psi_prev_val[i]) * uval[i] / TAU
    - H*H/(2*M) * (dPsi_iter_dx[i]*t_dudx + dPsi_iter_dy[i]*t_dudy)
    - G * Psi_iter_val[i] *  Psi_iter_val[i] * conj(Psi_iter_val[i]) * uval[i]
    - .5*M*OMEGA*OMEGA * (x[i]*x[i] + y[i]*y[i]) * Psi_iter_val[i] * uval[i]
  ));

  return result;
}

// Jacobian matrix for the implicit Euler time discretization
inline complex J_euler(ScalarFunction* Psi_iter, RealFunction* fu,
                RealFunction* fv, RefMap* ru, RefMap* rv)
{
  scalar ii = complex(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Quad2D* quad = fu->get_quad_2d();

  int o = 2 * Psi_iter->get_fn_order() + fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Psi_iter->set_quad_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  scalar* Psi_iter_val = Psi_iter->get_fn_values();
  scalar *dPsi_iter_dx, *dPsi_iter_dy;
  double* uval = fu->get_fn_values();
  double* vval = fv->get_fn_values();
  double *dudx, *dudy, *dvdx, *dvdy;
  Psi_iter->get_dx_dy_values(dPsi_iter_dx, dPsi_iter_dy);
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);

  // obtain physical coordinates of int. points
  double* x = ru->get_phys_x(o);
  double* y = ru->get_phys_y(o);

  // u is a basis function, v a test function
  scalar result = complex(0.0, 0.0);
  h1_integrate_dd_expression((
    ii * H * uval[i] * vval[i] / TAU
    - H*H/(2*M) * (t_dudx*t_dvdx + t_dudy*t_dvdy)
    - 2* G * uval[i] *  Psi_iter_val[i] * conj(Psi_iter_val[i]) * vval[i]
    - G * Psi_iter_val[i] *  Psi_iter_val[i] * uval[i] * vval[i]
    - .5*M*OMEGA*OMEGA * (x[i]*x[i] + y[i]*y[i]) * uval[i] * vval[i]
  ));

  return result;
}

// Residuum for the Crank-Nicolson time discretization
inline complex F_cranic(ScalarFunction* Psi_prev, ScalarFunction* Psi_iter, RealFunction* fu, RefMap* ru)
{
  scalar ii = complex(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Quad2D* quad = fu->get_quad_2d();
  RefMap* rv = ru;

  //not sure here:
  int o = 3 * Psi_iter->get_fn_order() + fu->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Psi_prev->set_quad_order(o);
  Psi_iter->set_quad_order(o);
  fu->set_quad_order(o);

  scalar* Psi_iter_val = Psi_iter->get_fn_values();
  scalar* Psi_prev_val = Psi_prev->get_fn_values();
  scalar* dPsi_iter_dx, *dPsi_iter_dy;
  scalar* dPsi_prev_dx, *dPsi_prev_dy;
  double* uval = fu->get_fn_values();
  double *dudx, *dudy;
  Psi_iter->get_dx_dy_values(dPsi_iter_dx, dPsi_iter_dy);
  Psi_prev->get_dx_dy_values(dPsi_prev_dx, dPsi_prev_dy);
  fu->get_dx_dy_values(dudx, dudy);

  // obtain physical coordinates of int. points
  double* x = ru->get_phys_x(o);
  double* y = ru->get_phys_y(o);

  // u is a test function
  scalar result = 0.0;
  h1_integrate_dd_expression((
    ii * H * (Psi_iter_val[i] - Psi_prev_val[i]) * uval[i] / TAU
    - 0.5*H*H/(2*M) * (dPsi_iter_dx[i]*t_dudx + dPsi_iter_dy[i]*t_dudy)
    - 0.5*H*H/(2*M) * (dPsi_prev_dx[i]*t_dudx + dPsi_prev_dy[i]*t_dudy)
    - 0.5*G * Psi_iter_val[i] *  Psi_iter_val[i] * conj(Psi_iter_val[i]) * uval[i]
    - 0.5*G * Psi_prev_val[i] *  Psi_prev_val[i] * conj(Psi_prev_val[i]) * uval[i]
    - 0.5*0.5*M*OMEGA*OMEGA * (x[i]*x[i] + y[i]*y[i]) * (Psi_iter_val[i] + Psi_prev_val[i]) * uval[i]
  ));

  return result;
}

// Jacobian matrix for the Crank-Nicolson time discretization
inline complex J_cranic(ScalarFunction* Psi_iter, RealFunction* fu,
                RealFunction* fv, RefMap* ru, RefMap* rv)
{
  scalar ii = complex(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Quad2D* quad = fu->get_quad_2d();

  int o = 2 * Psi_iter->get_fn_order() + fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  Psi_iter->set_quad_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  scalar* Psi_iter_val = Psi_iter->get_fn_values();
  scalar *dPsi_iter_dx, *dPsi_iter_dy;
  double* uval = fu->get_fn_values();
  double* vval = fv->get_fn_values();
  double *dudx, *dudy, *dvdx, *dvdy;
  Psi_iter->get_dx_dy_values(dPsi_iter_dx, dPsi_iter_dy);
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);

  // obtain physical coordinates of int. points
  double* x = ru->get_phys_x(o);
  double* y = ru->get_phys_y(o);

  // u is a basis function, v a test function
  scalar result = 0.0;
  h1_integrate_dd_expression((
    ii * H * uval[i] * vval[i] / TAU
    - 0.5*H*H/(2*M) * (t_dudx*t_dvdx + t_dudy*t_dvdy)
    - 0.5*2.0* G * uval[i] *  Psi_iter_val[i] * conj(Psi_iter_val[i]) * vval[i]
    - 0.5*G * Psi_iter_val[i] *  Psi_iter_val[i] * uval[i] * vval[i]
    - 0.5*0.5*M*OMEGA*OMEGA * (x[i]*x[i] + y[i]*y[i]) * uval[i] * vval[i]
  ));

  return result;
}
