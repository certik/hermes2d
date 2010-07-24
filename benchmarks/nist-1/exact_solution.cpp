static double fn(double x, double y)
{
  return pow(2.0,(4.0*p))*pow(x,p)*pow((1.0-x),p)*pow(y,p)*pow((1.0-y),p);
}

static double fndd(double x, double y, double& dx, double& dy)
{
  double A = pow((1.0-y),p);
  double B = pow((1.0-x),p);
  double C = pow(y,p);
  double D = pow(x,p);

  dx = ((p*pow(16.0,p)*A*C)/(x-1.0)+(p*pow(16.0,p)*A*C)/x)*B*D;
  dy = ((p*pow(16.0,p)*B*D)/(y-1.0)+(p*pow(16.0,p)*B*D)/y)*A*C;

  return fn(x, y);
}
