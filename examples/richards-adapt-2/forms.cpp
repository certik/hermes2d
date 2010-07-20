// Jacobian matrix - volumetric part
double jac_form_vol(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = ext->fn[0];
  Func<double>* h_prev_time = ext->fn[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * (
                         C(h_prev_newton->val[i]) * u->val[i] * v->val[i] / TAU
			 + dCdh(h_prev_newton->val[i]) * u->val[i] * h_prev_newton->val[i] * v->val[i] / TAU
                         - dCdh(h_prev_newton->val[i]) * u->val[i] * h_prev_time->val[i] * v->val[i] / TAU
			 + K(h_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                         + dKdh(h_prev_newton->val[i]) * u->val[i] * (h_prev_newton->dx[i]*v->dx[i] + h_prev_newton->dy[i]*v->dy[i])
                         - dKdh(h_prev_newton->val[i]) * u->dy[i] * v->val[i]
                         -  ddKdhh(h_prev_newton->val[i]) * u->val[i] * h_prev_newton->dy[i] * v->val[i]
                      );
  return result;
}

// Integration order for Jacobian matrix - volumetric part
Ord jac_form_vol_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Jacobian matrix - surface part on bdy 1
double jac_form_surf_1(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = ext->fn[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
  }

  return result;
}

// Integration order for Jacobian matrix - surface part on bdy 1
Ord jac_form_surf_1_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Jacobian matrix - surface part on bdy 4
double jac_form_surf_4(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = ext->fn[0];
  for (int i = 0; i < n; i++) {
    result -= wt[i] * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
  }

  return result;
}

// Integration order for Jacobian matrix - surface part on bdy 4
Ord jac_form_surf_4_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Jacobian matrix - surface part on bdy 6
double jac_form_surf_6(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = ext->fn[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
  }

  return result;
}

// Integration order for Jacobian matrix - surface part on bdy 6
Ord jac_form_surf_6_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Fesidual vector - volumetric part
double res_form_vol(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = ext->fn[0];
  Func<double>* h_prev_time = ext->fn[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * (
                          C(h_prev_newton->val[i]) * h_prev_newton->val[i] * v->val[i] / TAU
                        - C(h_prev_newton->val[i]) * h_prev_time->val[i] * v->val[i] / TAU
                        + K(h_prev_newton->val[i]) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
                        - dKdh(h_prev_newton->val[i]) * h_prev_newton->dy[i] * v->val[i]
                      );
  return result;
}


// Integration order for residual vector - volumetric part
Ord res_form_vol_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Fesidual vector - surface part on bdy 1
double res_form_surf_1(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = ext->fn[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * K(h_prev_newton->val[i]) * v->val[i];
  }

  return result;
}

// Integration order for residual vector - surface part on bdy 1
Ord res_form_surf_1_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Fesidual vector - surface part on bdy 4
double res_form_surf_4(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = ext->fn[0];
  for (int i = 0; i < n; i++) {
    result -= wt[i] * K(h_prev_newton->val[i]) * v->val[i];
  }

  return result;
}

// Integration order for residual vector - surface part on bdy 4
Ord res_form_surf_4_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Fesidual vector - surface part on bdy 6
double res_form_surf_6(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = ext->fn[0];
  for (int i = 0; i < n; i++) {
    result -= wt[i] * (Q_CONST - K(h_prev_newton->val[i])) * v->val[i];
  }

  return result;
}

// Integration order for residual vector - surface part on bdy 6
Ord res_form_surf_6_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

