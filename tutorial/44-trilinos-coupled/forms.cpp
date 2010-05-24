void omega_fn(int n, scalar* a, scalar* dadx, scalar* dady, scalar* b, scalar* dbdx, scalar* dbdy,
              scalar* out, scalar* outdx, scalar* outdy)
{
  for (int i = 0; i < n; i++)
  {
    double t1 = a[i] - 1.0;
    double t2 = t1 * beta;
    double t3 = 1.0 + t1 * alpha;
    double t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    double t5 = (beta / (t3 * t3)) * b[i];
    out[i] = t4 * b[i];
    outdx[i] = t4 * (dbdx[i] + dadx[i] * t5);
    outdy[i] = t4 * (dbdy[i] + dady[i] * t5);
  }
}

void omega_dt_fn(int n, scalar* a, scalar* dadx, scalar* dady, scalar* b, scalar* dbdx, scalar* dbdy,
                        scalar* out, scalar* outdx, scalar* outdy)
{
  for (int i = 0; i < n; i++)
  {
    double t1 = a[i] - 1.0;
    double t2 = t1 * beta;
    double t3 = 1.0 + t1 * alpha;
    double t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    double t5 = (beta / (t3 * t3));
    out[i] = t4 * t5 * b[i];
    outdx[i] = 0.0;
    outdy[i] = 0.0; // not important
  }
}

void omega_dc_fn(int n, scalar* a, scalar* dadx, scalar* dady, scalar* b, scalar* dbdx, scalar* dbdy,
                        scalar* out, scalar* outdx, scalar* outdy)
{
  for (int i = 0; i < n; i++)
  {
    double t1 = a[i] - 1.0;
    double t2 = t1 * beta;
    double t3 = 1.0 + t1 * alpha;
    double t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    out[i] = t4;
    outdx[i] = 0.0;
    outdy[i] = 0.0; // not important
  }
}

template<typename Real>
Real omega(Real t, Real c)
{
  Real t2 = beta * (t - 1.0);
  Real t3 = 1.0 + alpha * (t - 1.0);
  return c * sqr(beta) / (2.0*Le) * exp(t2 / t3);
}

template<typename Real, typename Scalar>
Scalar residual_0(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* tprev1 = ext->fn[0];
  Func<Real>* tprev2 = ext->fn[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( (3.0 * u[0]->val[i] - 4.0 * tprev1->val[i] + tprev2->val[i]) * vi->val[i] / (2.0 * TAU) +
                        (u[0]->dx[i] * vi->dx[i] + u[0]->dy[i] * vi->dy[i]) -
                        omega(u[0]->val[i], u[1]->val[i]) * vi->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar residual_0_surf(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (kappa * u[0]->val[i] * vi->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar residual_1(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* cprev1 = ext->fn[0];
  Func<Real>* cprev2 = ext->fn[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( (3.0 * u[1]->val[i] - 4.0 * cprev1->val[i] + cprev2->val[i]) * vi->val[i] / (2.0 * TAU) +
                        (u[1]->dx[i] * vi->dx[i] + u[1]->dy[i] * vi->dy[i]) / Le +
                        omega(u[0]->val[i], u[1]->val[i]) * vi->val[i] );
  return result;
}

// Preconditioner weak forms.
template<typename Real, typename Scalar>
Scalar precond_0_0(int n, double *wt, Func<Scalar>* u[], Func<Real> *vj, 
                   Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / TAU
                      +  vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar precond_1_1(int n, double *wt,  Func<Scalar>* u[], Func<Real> *vj, 
                   Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / TAU
                      +  (vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]) / Le );
  return result;
}

// Jacobian weak forms.
template<typename Real, typename Scalar>
Scalar jacobian_0_0(int n, double *wt, Func<Scalar>* u[], Func<Real> *vj, 
                    Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / TAU
                      +  vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]
                      - (u[1]->val[i] * sqr(beta)/(2*Le) *
                         exp(beta * (u[0]->val[i] - 1)/(1.0 + alpha*(u[0]->val[i] - 1))) *
                         beta / (sqr(1.0 + alpha*(u[0]->val[i] - 1))))
                           * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar jacobian_0_0_surf(int n, double *wt, Func<Scalar>* u[], Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (kappa * vj->val[i] * vi->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar jacobian_0_1(int n, double *wt, Func<Scalar>* u[], Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* dodc = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (- (sqr(beta)/(2*Le) * exp(beta * (u[0]->val[i]-1)/(1.0 + alpha*(u[0]->val[i]-1)))) * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar jacobian_1_0(int n, double *wt, Func<Scalar>* u[], Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* dodt = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( (u[1]->val[i] * sqr(beta)/(2*Le) *
                         exp(beta * (u[0]->val[i] - 1)/(1.0 + alpha*(u[0]->val[i] - 1))) *
                         beta / (sqr(1.0 + alpha*(u[0]->val[i] - 1)))) * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar jacobian_1_1(int n, double *wt,  Func<Scalar>* u[], Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / TAU
                      +  (vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]) / Le
                      + (sqr(beta)/(2*Le) * exp(beta * (u[0]->val[i]-1)/(1.0 + alpha*(u[0]->val[i]-1))))
                        * vj->val[i] * vi->val[i] );
  return result;
}
