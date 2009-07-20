#include "hermes2d.h"
#include "solver_umfpack.h"

//
//  PDE: time-harmonic Maxwell's equations
//
//  Known exact solution, see the function exact()
//
//  Domain: square domain cut from the midpoint of the left edge to the center (center
//          point of left edge duplicated)
//
//  Meshes: you can either use "screen-quad.mesh" (quadrilateral mesh) or
//          "screen-tri.mesh" (triangular mesh). See the command mesh.load(...) below
//
//  BC: tangential component of solution taken from known exact solution (essential BC),
//      see function bc_values(...) below
//
//  This example has a known exact solution. It describes an electromagnetic wave that hits
//  a screen under the angle of 45 degrees, causing a singularity at the tip of the screen.
//  Convergence graphs saved (both exact error and error estimate, and both wrt. dof number
//  and cpu time). The following problem parameters can be changed easily:
//

int P_INIT = 1;           // initial polynomial degree in mesh
double THRESHOLD = 0.5;   // the adaptivity algorithm goes on until THRESHOLD*total_error is processed
                          // (see adapt_hcurl.cpp for explanation)
int STRATEGY = 1;         // refinement strategy (0, 1, 2, 3 - see adapt_hcurl.cpp for explanation)
int H_ONLY = 0;           // if H_ONLY == 0 then full hp-adaptivity takes place, otherwise
                          // h-adaptivity is used. Use this parameter to check that indeed adaptive
                          // hp-FEM converges much faster than adaptive h-FEM
double ERR_STOP = 0.5;    // adaptivity process stops when error wrt. error estimate in H(curl) norm
                          // is less than this number
int NDOF_STOP = 40000;    // adaptivity process stops when the number of degrees of freedom grows over
                          // this limit. This is mainly to prevent h-adaptivity to go on forever.

const double e_0  = 8.8541878176 * 1e-12;
const double mu_0 = 1.256 * 1e-6;
const double k = 1.0;


//// exact solution ////////////////////////////////////////////////////////////////////////////////

extern "C" void fresnl( double xxa, double *ssa, double *cca );


scalar Fn(double u)
{
  double s, c;
  fresnl(sqrt(2/M_PI) * u, &s , &c);
  scalar fres = complex(c,-s);
  scalar a = complex(0.0, M_PI/4);
  scalar b = complex(0.0, u*u);
  return 0.5*sqrt(M_PI) * exp(b) * (exp(-a) - sqrt(2)*(fres));
}

scalar Fder(double u)
{
  scalar a = complex(0.0, M_PI/4);
  scalar b = complex(0.0, u*u);
  scalar d = complex(0.0, 2.0*u);
  double s, c;
  fresnl(sqrt(2/M_PI) * u, &s , &c);
  scalar fres = complex(c,-s);
  scalar fresder = exp(-b);

  return 0.5*sqrt(M_PI) * exp(b) * ( d * (exp(-a) - sqrt(2)*(fres)) - sqrt(2)*fresder*sqrt(2.0/M_PI) );
}

scalar Fder2(double u)
{
  scalar a = complex(0.0, M_PI/4);
  scalar i = complex(0.0,1.0);
  scalar b = complex(0.0, u*u);
  scalar d = complex(0.0, 2.0*u);
  double s, c;
  fresnl(sqrt(2/M_PI) * u, &s , &c);
  scalar fres = complex(c,-s);
  scalar fresder = exp(-b);
  scalar fresder2 = exp(-b)*(-2.0 * i * u);

  return 2.0 * u * i * Fder(u) +
         0.5 * sqrt(M_PI) * exp(b) *
          ( 2.0 * i * (exp(-a) - sqrt(2)*(fres)) + d * (-sqrt(2)*fresder*sqrt(2.0/M_PI)) - sqrt(2) * fresder2 * sqrt(2.0/M_PI) );
}

scalar der_Hr(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar a = complex(0.0, M_PI/4 - k*r);
  scalar i = complex(0.0,1.0);
  return 1/sqrt(M_PI) * exp(a) *
        ( (-i*k)*(Fn(sqrt(2*k*r)*sin(t/2 - M_PI/8)) + Fn(sqrt(2*k*r)*sin(t/2 + M_PI/8))) +
        (Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8))*(sqrt(k)/sqrt(2*r)*sin(t/2 - M_PI/8)) +
         Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8))*(sqrt(k)/sqrt(2*r)*sin(t/2 + M_PI/8))));
}

scalar der_Hrr(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar a = complex(0.0, M_PI/4 - k*r);
  scalar i = complex(0.0,1.0);
  scalar f1_d = Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d = Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar f1_d2 = Fder2(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d2 = Fder2(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar b1 = (sqrt(k/(2*r))*sin(t/2 - M_PI/8));
  scalar b2 = (sqrt(k/(2*r))*sin(t/2 + M_PI/8));
  return -i*k*der_Hr(x,y) + 1/sqrt(M_PI) * exp(a) *
        ( (-i*k)*(f1_d*b1 + f2_d*b2) +
        ( f1_d2*b1*b1 + f2_d2*b2*b2) +
          f1_d*(-0.5*sqrt(k/(2*r*r*r))*sin(t/2 - M_PI/8))  + f2_d*(-0.5*sqrt(k/(2*r*r*r))*sin(t/2 + M_PI/8)));
}

scalar der_Hrt(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar i = complex(0.0,1.0);
  scalar a = complex(0.0, M_PI/4 - k*r);
  scalar f1_d = Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d = Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar f1_d2 = Fder2(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d2 = Fder2(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar b1 = (sqrt(k)/sqrt(2*r)*sin(t/2 - M_PI/8));
  scalar b2 = (sqrt(k)/sqrt(2*r)*sin(t/2 + M_PI/8));
  scalar c1 = (sqrt(k*r)/sqrt(2)*cos(t/2 - M_PI/8));
  scalar c2 = (sqrt(k*r)/sqrt(2)*cos(t/2 + M_PI/8));
  return 1/sqrt(M_PI) * exp(a) *
        ( (-i*k)*(f1_d*c1 + f2_d*c2) +
        ( f1_d2*b1*c1 + f2_d2*b2*c2) +
          f1_d*(0.5*sqrt(k/(2*r))*cos(t/2 - M_PI/8))  + f2_d*(0.5*sqrt(k/(2*r))*cos(t/2 + M_PI/8)));
}

scalar der_Ht(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar a = complex(0.0, M_PI/4 - k*r);
  return 1/sqrt(M_PI) * exp(a) *
         (Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8))*(sqrt(k*r/2)*cos(t/2 - M_PI/8)) +
          Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8))*(sqrt(k*r/2)*cos(t/2 + M_PI/8)));
}

scalar der_Htr(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar i = complex(0.0,1.0);
  scalar a = complex(0.0, M_PI/4 - k*r);
  scalar f1_d = Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d = Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar f1_d2 = Fder2(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d2 = Fder2(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar b1 = (sqrt(k)/sqrt(2*r)*sin(t/2 - M_PI/8));
  scalar b2 = (sqrt(k)/sqrt(2*r)*sin(t/2 + M_PI/8));
  scalar c1 = (sqrt(k*r)/sqrt(2)*cos(t/2 - M_PI/8));
  scalar c2 = (sqrt(k*r)/sqrt(2)*cos(t/2 + M_PI/8));
  return -i*k*der_Ht(x,y) + 1/sqrt(M_PI) * exp(a) *
         ((f1_d2*b1*c1 + f2_d2*b2*c2) +
          f1_d*(0.5*sqrt(k/(2*r))*cos(t/2 - M_PI/8))  + f2_d*(0.5*sqrt(k/(2*r))*cos(t/2 + M_PI/8)));
}

scalar der_Htt(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar a = complex(0.0, M_PI/4 - k*r);
  scalar f1_d = Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d = Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar f1_d2 = Fder2(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d2 = Fder2(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar c1 = (sqrt(k*r/(2))*cos(t/2 - M_PI/8));
  scalar c2 = (sqrt(k*r/(2))*cos(t/2 + M_PI/8));
  return 1/sqrt(M_PI) * exp(a) *
         ((f1_d2*c1*c1 + f2_d2*c2*c2) +
          f1_d*(-0.5*sqrt(k*r/2)*sin(t/2 - M_PI/8))  + f2_d*(-0.5*sqrt(k*r/2)*sin(t/2 + M_PI/8)));
}

scalar exact0(double x, double y, scalar& dx, scalar& dy)
{
  double r = sqrt(x*x + y*y);
  double theta = atan2(y,x);
  scalar Hr = der_Hr(x,y);
  scalar Ht = der_Ht(x,y);
  scalar i = complex(0.0,1.0);
  return  -i * (Hr * y/r + Ht * x/(r*r));
}

scalar exact1(double x, double y, scalar& dx, scalar& dy)
{
  double r = sqrt(x*x + y*y);
  double theta = atan2(y,x);
  scalar Hr = der_Hr(x,y);
  scalar Ht = der_Ht(x,y);
  scalar i = complex(0.0,1.0);
  return  i * ( Hr * x/r - Ht * y/(r*r));
}

static void exact_sol(double x, double y, scalar& u0, scalar& u1, scalar& u1dx, scalar& u0dy)
{
  scalar dx,dy;
  u0 = exact0(x,y,dx,dy);
  u1 = exact1(x,y,dx,dy);

  scalar Hr = der_Hr(x,y);
  scalar Ht = der_Ht(x,y);
  scalar Hrr = der_Hrr(x,y);
  scalar Hrt = der_Hrt(x,y);
  scalar Htr = der_Htr(x,y);
  scalar Htt = der_Htt(x,y);

  double r = sqrt(x*x + y*y);
  double theta = atan2(y,x);
  scalar i = complex(0.0,1.0);

  u1dx =  i * (( Hrr * x/r + Hrt * (-y/(r*r))) * x/r     + Hr * (y*y)/(r*r*r) -
               ((Htr * x/r + Htt * (-y/(r*r))) * y/(r*r) + Ht * (-2.0*x*y/(r*r*r*r))));
  u0dy = -i * (( Hrr * y/r + Hrt *   x/(r*r))  * y/r     + Hr * (x*x)/(r*r*r) +
                (Htr * y/r + Htt *   x/(r*r))  * x/(r*r) + Ht * (-2.0*x*y/(r*r*r*r)));
}

scalar2& exact(double x, double y, scalar2& dx, scalar2& dy)
{
  scalar2 ex;
  exact_sol(x,y, ex[0], ex[1], dx[1], dy[0]);

  dx[0] = 0.0; // not important
  dy[1] = 0.0; // not important

  return ex;
}


//// problem definition ////////////////////////////////////////////////////////////////////////////////////////////

int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// TODO: obtain tangent from EdgePos
double2 tau[5] = { { 0, 0}, { 1, 0 },  { 0, 1 }, { -1, 0 }, { 0, -1 } };

complex bc_values(int marker, double x, double y)
{
  scalar dx, dy;
  return exact0(x, y, dx, dy)*tau[marker][0] + exact1(x, y, dx, dy)*tau[marker][1];
}

complex bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_curl_e_curl_f(fu, fv, ru, rv) - int_e_f(fu, fv, ru, rv);
}


//// main //////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("screen-quad.mesh");
  //mesh.load("screen-tri.mesh");

  HcurlShapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  HcurlSpace space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form, SYM);

  ScalarView Xview_r("Electric field X - real",   0, 0, 320, 320);
  ScalarView Yview_r("Electric field Y - real", 325, 0, 320, 320);
  ScalarView Xview_i("Electric field X - imag", 650, 0, 320, 320);
  ScalarView Yview_i("Electric field Y - imag", 975, 0, 320, 320);

  OrderView  ord("Polynomial Orders", 325, 400, 600, 600);

  GnuplotGraph graph;
  graph.set_captions("Error Convergence for the Screen Problem in H(curl)", "Degrees of Freedom", "Error [%]");
  graph.add_row("exact error", "k", "-", "o");
  graph.add_row("error estimate", "k", "--");
  graph.set_log_y();

  GnuplotGraph graph_cpu;
  graph_cpu.set_captions("Error Convergence for the Screen Problem in H(curl)", "CPU Time", "Error [%]");
  graph_cpu.add_row("exact error", "k", "-", "o");
  graph_cpu.add_row("error estimate", "k", "--");
  graph_cpu.set_log_y();

  UmfpackSolver umfpack;
  Solution sln, rsln;

  int it = 0;
  begin_time();
  bool done = false;
  double cpu = 0.0;
  do
  {
    printf("\n---- it=%d ------------------------------------------------------------------\n\n", it++);

    // enumerating basis functions
    space.assign_dofs();
    begin_time();

    // coarse problem
    LinSystem sys(&wf, &umfpack);
    sys.set_spaces(1, &space);
    sys.set_pss(1, &pss);
    sys.assemble();
    sys.solve(1, &sln);

    cpu += end_time();

    // visualization
    RealFilter real(&sln);
    ImagFilter imag(&sln);
    Xview_r.set_min_max_range(-3.0, 1.0);
    Xview_r.show_scale(false);
    Xview_r.show(&real, EPS_NORMAL, FN_VAL_0);
    Yview_r.set_min_max_range(-4.0, 4.0);
    Yview_r.show_scale(false);
    Yview_r.show(&real, EPS_NORMAL, FN_VAL_1);
    Xview_i.set_min_max_range(-1.0, 4.0);
    Xview_i.show_scale(false);
    Xview_i.show(&imag, EPS_NORMAL, FN_VAL_0);
    Yview_i.set_min_max_range(-4.0, 4.0);
    Yview_i.show_scale(false);
    Yview_i.show(&imag, EPS_NORMAL, FN_VAL_1);

    ord.show(&space);

    // calculating error wrt. exact solution
    ExactSolution ex(&mesh, exact);
    double error = 100 * hcurl_error(&sln, &ex);

    // fine (reference) problem
    begin_time();
    RefSystem ref(&sys);
    ref.assemble();
    ref.solve(1, &rsln);

    // calculating error estimate wrt. fine mesh solution
    HcurlOrthoHP hp(1, &space);
    double estim = hp.calc_error(&sln, &rsln) * 100;
    info("Exact solution error: %g%%", error);
    info("Error estimate: %g%%", estim);

    // mesh adaptation
    if (estim < ERR_STOP || sys.get_num_dofs() >= NDOF_STOP) done = true;
    hp.adapt(THRESHOLD, STRATEGY, H_ONLY);

    // plotting convergence wrt. numer of dofs
    graph.add_values(0, space.get_num_dofs(), error);
    graph.add_values(1, space.get_num_dofs(), estim);
    graph.save("conv_dof.gp");

    // plotting convergence wrt. numer of dofs
    cpu += end_time();
    graph_cpu.add_values(0, cpu, error);
    graph_cpu.add_values(1, cpu, estim);
    graph_cpu.save("conv_cpu.gp");
  }
  while (!done);
  verbose("\nTotal run time: %g sec", end_time());

  View::wait();
  return 0;
}

