static double fn(double x, double y)
{
  return atan(ALPHA * (sqrt(pow(x - X_LOC, 2) + pow(y - Y_LOC, 2)) - R_ZERO));
}

static double fndd(double x, double y, double& dx, double& dy)
{
  double a = pow(x - X_LOC, 2);
  double b = pow(y - Y_LOC, 2);
  double c = sqrt(a + b);
  double d = (ALPHA*x - 25);
  double e = (ALPHA*y - 25);
  double f = (pow(ALPHA * c - 12.5, 2) + 1);

  dx = (d/(f * c));
  dy = (e/(f * c));

  return fn(x, y);
}
