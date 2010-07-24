template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real>
Real rhs(Real x, Real y)
{
  Real A = pow(2.0,(4.0*p));
  Real B = pow((x-1.0),8.0);
  Real C = (38.0*pow(x,2.0)-38.0*x+9.0);
  Real D = pow((y-1.0),p);
  Real E = pow((y-1.0),8.0);
  Real F = (38.0*pow(y,2.0)-38.0*y+9.0);
  Real G = pow((x-1.0),p);

  return p*A*pow(x,8.0)*B*C*pow(y,p)*D + p*A*pow(y,8.0)*E*F*pow(x,p)*G;
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return -int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}
