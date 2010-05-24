template<typename Real, typename Scalar>
Scalar jacobian_form_hermes(int n, double *wt, Func<Real> *vi, Func<Real> *vj, 
                            Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* u = ext->fn[0];
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( -0.5 * pow(1.0 + sqr(u->dx[i]) + sqr(u->dy[i]), -1.5) * 
                       (2.0 * u->dx[i] * vi->dx[i] + 2.0 * u->dx[i] * vi->dx[i])
                       * (u->dx[i] * vj->dx[i] + u->dy[i] * vj->dy[i]) +
                       (pow(1.0 + sqr(u->dx[i]) + sqr(u->dy[i]), -0.5))
                       * (vi->dx[i] * vj->dx[i] + vi->dy[i] * vj->dy[i]) );
  return result;
}

template<typename Real, typename Scalar>
Scalar residual_form_hermes(int n, double *wt, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* u = ext->fn[0];
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ((pow(1.0 + sqr(u->dx[i]) + sqr(u->dy[i]), -0.5)) * (u->dx[i] * vj->dx[i] + u->dy[i] * vj->dy[i])
                       - f(e->x[i], e->y[i]) * vj->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar jacobian_form_nox(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( -0.5 * pow(1.0 + sqr(u[0]->dx[i]) + sqr(u[0]->dy[i]), -1.5) *
                              (2.0 * u[0]->dx[i] * vi->dx[i] + 2.0 * u[0]->dx[i] * vi->dx[i])
                       * (u[0]->dx[i] * vj->dx[i] + u[0]->dy[i] * vj->dy[i]) +
                       (pow(1.0 + sqr(u[0]->dx[i]) + sqr(u[0]->dy[i]), -0.5))
                       * (vi->dx[i] * vj->dx[i] + vi->dy[i] * vj->dy[i]) );
  return result;
}

template<typename Real, typename Scalar>
Scalar precond_form_nox(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( vi->dx[i] * vj->dx[i] + vi->dy[i] * vj->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar residual_form_nox(int n, double *wt, Func<Real> *u[], Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ((pow(1.0 + sqr(u[0]->dx[i]) + sqr(u[0]->dy[i]), -0.5)) *
                        (u[0]->dx[i] * vj->dx[i] + u[0]->dy[i] * vj->dy[i])
                       - f(e->x[i], e->y[i]) * vj->val[i] );
  return result;
}
