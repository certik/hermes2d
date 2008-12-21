#include "hermes2d.h"
#include "solver_umfpack.h"

const double mu_r   = 1.0;
const double kappa  = 1.0;
const double lambda = 1.0;

// Bessel function of the first kind, order n, defined in bessel.cpp
double jv(double n, double x);


//// exact solution ////////////////////////////////////////////////////////////////////////////////

static void exact_sol_val(double x, double y, scalar& e0, scalar& e1)
{
  double t1 = x*x;
  double t2 = y*y;
  double t4 = sqrt(t1+t2);
  double t5 = jv(-1.0/3.0,t4);
  double t6 = 1/t4;
  double t7 = jv(2.0/3.0,t4);
  double t11 = (t5-2.0/3.0*t6*t7)*t6;
  double t12 = atan2(y,x);
  if (t12 < 0) t12 += 2.0*M_PI;
  double t13 = 2.0/3.0*t12;
  double t14 = cos(t13);
  double t17 = sin(t13);
  double t18 = t7*t17;
  double t20 = 1/t1;
  double t23 = 1/(1.0+t2*t20);
  e0 = t11*y*t14-2.0/3.0*t18/x*t23;
  e1 = -t11*x*t14-2.0/3.0*t18*y*t20*t23;
}

static void exact_sol(double x, double y, scalar& e0, scalar& e1, scalar& e1dx, scalar& e0dy)
{
  exact_sol_val(x,y,e0,e1);
  
  double t1 = x*x;
  double t2 = y*y;
  double t3 = t1+t2;
  double t4 = sqrt(t3);
  double t5 = jv(2.0/3.0,t4);
  double t6 = 1/t4;
  double t7 = jv(-1.0/3.0,t4);
  double t11 = (-t5-t6*t7/3.0)*t6;
  double t14 = 1/t4/t3;
  double t15 = t14*t5;
  double t21 = t7-2.0/3.0*t6*t5;
  double t22 = 1/t3*t21;
  double t27 = atan2(y,x);
  if (t27 < 0) t27 += 2.0*M_PI;
  double t28 = 2.0/3.0*t27;
  double t29 = cos(t28);
  double t32 = t21*t14;
  double t35 = t21*t6;
  double t36 = t35*t29;
  double t39 = sin(t28);
  double t41 = 1/t1;
  double t43 = 1.0+t2*t41;
  double t44 = 1/t43;
  double t47 = 4.0/3.0*t35/x*t39*y*t44;
  double t48 = t5*t29;
  double t49 = t1*t1;
  double t52 = t43*t43;
  double t53 = 1/t52;
  double t57 = t5*t39;
  double t59 = 1/t1/x;
  e1dx =-(t11*x+2.0/3.0*t15*x-2.0/3.0*t22*x)
              *t6*x*t29+t32*t1*t29-t36-t47+4.0/9.0*t48*t2/t49*t53+4.0/3.0*t57*y*t59*t44-4.0/3.0*t57*t2*y/t49/x*t53;
  e0dy = (t11*y+2.0/3.0*t15*y-2.0/3.0*t22*y)*t6*y*t29-t32*t2*t29+t36-t47-4.0/9.0*t48*t41*t53+4.0/3.0*t57*t59*t53*y;
}


scalar2& exact(double x, double y, scalar2& dx, scalar2& dy)
{
  scalar2 ex;
  exact_sol(x,y, ex[0], ex[1], dx[1], dy[0]);
  return ex;
}


//// boundary condition ////////////////////////////////////////////////////////////////////////////

int bc_types(int marker)
{
  if (marker == 1 || marker == 6)
    return BC_ESSENTIAL; // perfect conductor
  else
    return BC_NATURAL; // impedance
}

// TODO: obtain tangent from EdgePos
double2 tau[7] = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 }, { 1, 0 }, { 0, 1 } };
  
complex bc_values(int marker, double x, double y)
{
  if (marker == 1 || marker == 6) return 0;

  double r = sqrt(x*x + y*y), theta = atan2(y, x);
  if (theta < 0) theta += 2.0*M_PI;
  double j13    = jv(-1.0/3.0, r),    j23    = jv(+2.0/3.0, r);
  double cost   = cos(theta),         sint   = sin(theta);
  double cos23t = cos(2.0/3.0*theta), sin23t = sin(2.0/3.0*theta);

  double Etau = tau[marker][0] * (cos23t*sint*j13 - 2.0/(3.0*r)*j23*(cos23t*sint + sin23t*cost)) +
                tau[marker][1] * (-cos23t*cost*j13 + 2.0/(3.0*r)*j23*(cos23t*cost - sin23t*sint));

  return complex(cos23t*j23, -Etau);
}


//// weak formulation //////////////////////////////////////////////////////////////////////////////

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


//// main //////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("lshape3q.mesh");
  //mesh.load("lshape3t.mesh");

  HcurlShapeset shapeset;
  PrecalcShapeset pss(&shapeset);
 
  HcurlSpace space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(1);

  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form, SYM);
  wf.add_biform_surf(0, 0, bilinear_form_surf);
  wf.add_liform_surf(0, linear_form_surf);
  
  GnuplotGraph graph;
  graph.set_captions("Error Convergence for the Bessel Problem in H(curl)", "Degrees of Freedom", "Error [%]");
  graph.add_row("ortho adaptivity", "k", "-", "o");
  graph.set_log_y();

  OrderView  ord("Polynomial Orders", 800, 100, 700, 600); 
  VectorView vecview("Real part of Electric Field - VectorView", 0, 100, 700, 600);

  UmfpackSolver umfpack;
  Solution sln, rsln;

  int it = 0;
  begin_time();
  while (1)
  {
    printf("\n\n---- it=%d -----------------------------------------------------------\n\n", it++);

    space.assign_dofs();
    ord.show(&space);
    
    // coarse problem
    LinSystem sys(&wf, &umfpack);
    sys.set_spaces(1, &space);
    sys.set_pss(1, &pss);
    sys.assemble();
    sys.solve(1, &sln);

    // get the exact error of the solution
    ExactSolution ex(&mesh, exact);
    double error = 100 * hcurl_error(&sln, &ex);
    info("\nExact solution error: %g%%\n", error);
    graph.add_values(0, space.get_num_dofs(), error);
    graph.save("convergence.gp");
    
    // show real part of the solution
    RealFilter real(&sln);
    vecview.set_min_max_range(0, 1);
    vecview.show(&real, EPS_HIGH);
    
    // fine (reference) problem
    RefSystem ref(&sys);
    ref.assemble();
    ref.solve(1, &rsln);

    HcurlOrthoHP hp(1, &space);
    double estim = hp.calc_error(&sln, &rsln) * 100;
    info("\nError estimate: %g%%", estim);
    if (estim < 0.01) break;
    
    // adapt the mesh&space and repeat
    hp.adapt(0.3, 1);
  }
  verbose("\nTotal running time: %g sec", end_time());

  View::wait();
  return 0;
}

