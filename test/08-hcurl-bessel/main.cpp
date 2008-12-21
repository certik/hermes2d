#include "hermes2d.h"

const double mu_r   = 1.0;
const double kappa  = 1.0;
const double lambda = 1.0;

// Bessel function of the first kind, order n, defined in bessel.cpp
double jv(double n, double x);


int bc_types(int marker)
{
  if (marker == 1 || marker == 6)
    return BC_ESSENTIAL; // perfect conductor
  else
    return BC_NATURAL; // impedance
}


double2 tau[7] = { { 0, 0 }, { -1, 0 }, { 0, -1 }, { -1, 0 }, { 0, 1 }, { 1, 0 }, { 0, -1 } };
  
complex bc_values(int marker, double x, double y)
{
  // perfect conductor BC: E.tau = 0
  if (marker == 1 || marker == 6) return 0;

  // impedance BC: return g.tau
  double r = sqrt(x*x + y*y), theta = atan2(y, x);
  if (theta < 0) theta += 2.0*M_PI;
  double j13    = jv(-1.0/3.0, r),    j23    = jv(+2.0/3.0, r);
  double cost   = cos(theta),         sint   = sin(theta);
  double cos23t = cos(2.0/3.0*theta), sin23t = sin(2.0/3.0*theta);

  double Etau = tau[marker][0] * (cos23t*sint*j13 - 2.0/(3.0*r)*j23*(cos23t*sint + sin23t*cost)) +
                tau[marker][1] * (cos23t*cost*j13 + 2.0/(3.0*r)*j23*(cos23t*cost - sin23t*sint));

  return complex(cos23t*j23, -Etau);
}


complex bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return    1.0/mu_r * int_curl_e_curl_f(fu, fv, ru, rv)
        - sqr(kappa) * int_e_f(fu, fv, ru, rv);
}

complex bilinear_form_surf(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv, EdgePos* ep)
{
  if (ep->marker == 1 || ep->marker == 6) return 0;
  return complex(0.0, -kappa * surf_int_e_tau_f_tau(fu, fv, ru, rv, ep));
}

complex linear_form_surf(RealFunction* fv, RefMap* refmap, EdgePos* ep)
{
  if (ep->marker == 1 || ep->marker == 6) return 0;
  return surf_int_G_tau_f_tau(fv, refmap, ep);
}


scalar exact0(double x, double y, scalar& dx, scalar& dy)
{
  double r = sqrt(x*x + y*y), theta = atan2(y, x);
  if (theta < 0) theta += 2.0*M_PI;
  double j13    = jv(-1.0/3.0, r),    j23    = jv(+2.0/3.0, r);
  double cost   = cos(theta),         sint   = sin(theta);
  double cos23t = cos(2.0/3.0*theta), sin23t = sin(2.0/3.0*theta);
  return cos23t*sint*j13 - 2.0/(3.0*r)*j23*(cos23t*sint + sin23t*cost);
}

scalar exact1(double x, double y, scalar& dx, scalar& dy)
{
  double r = sqrt(x*x + y*y), theta = atan2(y, x);
  if (theta < 0) theta += 2.0*M_PI;
  double j13    = jv(-1.0/3.0, r),    j23    = jv(+2.0/3.0, r);
  double cost   = cos(theta),         sint   = sin(theta);
  double cos23t = cos(2.0/3.0*theta), sin23t = sin(2.0/3.0*theta);
  return -cos23t*cost*j13 + 2.0/(3.0*r)*j23*(cos23t*cost - sin23t*sint);
}


int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);
  
  Mesh mesh;
  mesh.load("lshape3q.mesh");
  mesh.refine_towards_vertex(0, 5);
  mesh.refine_all_elements();
  
  HcurlShapesetGradLeg shapeset;
  PrecalcShapeset pss(&shapeset);

  HcurlSpace space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(3);
  space.assign_dofs();

  DiscreteProblem dp;
  dp.set_num_equations(1);
  dp.set_spaces(1, &space);
  dp.set_pss(1, &pss);
  dp.set_bilinear_form(0, 0, bilinear_form, NULL, bilinear_form_surf);
  dp.set_linear_form(0, NULL, linear_form_surf);

  Solution sln;
  dp.create_matrix();
  dp.assemble_matrix_and_rhs();
  dp.solve_system(1, &sln);
  
  sln.save("bessel.sln");
 
  RealFilter real(&sln);
  ScalarView view1("X component", 100, 150, 1000, 900);
  view1.show(&real, EPS_NORMAL, FN_VAL_0);
  ScalarView view2("Y component", 200, 150, 1000, 900);
  view2.show(&real, EPS_NORMAL, FN_VAL_1);

  MagFilter mag(&real);
  ScalarView view3("Magnitude of real(E)", 100, 50, 1000, 900);
  view3.set_min_max_range(0, 1);
  view3.show_contours(0.07);
  view3.show(&mag, EPS_HIGH);

  hermes2d_finalize();
  return 0;
}
