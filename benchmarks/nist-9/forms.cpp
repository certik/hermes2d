template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real>
Real rhs(Real x, Real y)
{
  Real a = pow(x - X_LOC, 2);
  Real b = pow(y - Y_LOC, 2);
  Real c = sqrt(a + b);
  Real d = ((2 * x - 1.0) * (ALPHA * x - 25.0));
  Real e = ((2 * y - 1.0) * (ALPHA * y - 25.0));
  Real f = (pow(ALPHA * c - 12.5, 2) + 1);

  return ((ALPHA/(f * c)) - (d/(2 * f * pow(a + b, 1.5))) - ((ALPHA * (ALPHA * c - 12.5) * d)/ (pow(f, 2) * (a + b))) +
         (ALPHA/(f * c)) - (e/(2 * f * pow(a + b, 1.5))) - ((ALPHA * (ALPHA * c - 12.5) * e)/ (pow(f, 2) * (a + b))));
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return -int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}
