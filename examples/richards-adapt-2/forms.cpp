// Jacobian matrix.
double jac(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext)
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

Ord jac_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Fesidual vector
double res(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
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

Ord res_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}
