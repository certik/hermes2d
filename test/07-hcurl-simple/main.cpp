#include "hermes2d.h"

const double mu_r   = 1.0;
const double kappa  = 1.0;
const double lambda = 1.0;


int bc_types(int marker)
{
  if (marker == 1 || marker == 6)
    return BC_ESSENTIAL; // perfect conductor
  else
    return BC_NATURAL; // impedance
}

double2 tau[7] = { { 0, 0 }, { 0, -1 }, { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 }, { 1, 0 } };

complex bc_values(int marker, double x, double y)
{
  // perfect conductor BC: E.tau = 0
  if (marker == 1 || marker == 6) return 0;

  // impedance BC: return g.tau
  double theta = atan2(y, x);
  double r13 = pow(sqrt(x*x + y*y), -1.0/3.0);
  return complex(0.0, tau[marker][0] * (-2.0/3.0) * r13 * cos(M_PI/6 + theta/3) +
                      tau[marker][1] * (-2.0/3.0) * r13 * sin(M_PI/6 + theta/3));
}

void rhs(double x, double y, complex2& result)
{
  double theta = atan2(y, x);
  double r13 = pow(sqrt(x*x + y*y), -1.0/3.0);
  result[0] = (-2.0/3.0) * r13 * cos(M_PI/6 + theta/3);
  result[1] = (-2.0/3.0) * r13 * sin(M_PI/6 + theta/3);
}

complex bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return    1.0/mu_r * int_curl_e_curl_f(fu, fv, ru, rv)
        - sqr(kappa) * int_e_f(fu, fv, ru, rv);
}

complex bilinear_form_surf(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv, EdgePos* ep)
{
  if (ep->marker == 1 || ep->marker == 6) return 0;
  return complex(0.0, -surf_int_e_tau_f_tau(fu, fv, ru, rv, ep));
}

complex linear_form(RealFunction* fv, RefMap* rv)
{
  return int_F_f(rhs, fv, rv);
}

complex linear_form_surf(RealFunction* fv, RefMap* rv, EdgePos* ep)
{
  if (ep->marker == 1 || ep->marker == 6) return 0;
  return surf_int_G_tau_f_tau(fv, rv, ep);
}


int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);

  Mesh mesh;
  mesh.load("lshape2.mesh");
  mesh.refine_towards_vertex(0, 5);

  HcurlShapesetLegendre shapeset;
  PrecalcShapeset pss(&shapeset);

  HcurlSpace space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(8);
  space.assign_dofs();

  DiscreteProblem dp;
  dp.set_num_equations(1);
  dp.set_spaces(1, &space);
  dp.set_pss(1, &pss);
  dp.set_bilinear_form(0, 0, NULL, bilinear_form, bilinear_form_surf);
  dp.set_linear_form(0, linear_form, linear_form_surf);

  Solution sln;
  dp.create_matrix();
  dp.assemble_matrix_and_rhs();
  dp.solve_system(1, &sln);
  
  RealFilter real(&sln);
  
  ScalarView view1("X component", 100, 50, 900, 900);
  view1.show(&real, EPS_HIGH, FN_VAL_0);
  view1.set_min_max_range(0, 5);
  
  ScalarView view2("Y component", 200, 150, 900, 900);
  view2.show(&real, EPS_HIGH, FN_VAL_1);
  view2.set_min_max_range(0, 5);

  MagFilter mag(&real, &real, FN_VAL_0, FN_VAL_1);
  ScalarView view3("Magnitude of E", 300, 250, 900, 900);
  view3.set_min_max_range(0.2, 2.5);
  view3.show_contours(0.1);
  view3.show(&mag, EPS_HIGH);
  
  VectorView view4;
  view4.show(&sln);

  hermes2d_finalize();
  return 0;
}
