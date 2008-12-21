#include "hermes2d.h"
#include "solver_umfpack.h"


static double fn(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double a = atan2(x, y);
  return pow(r, 2.0/3.0) * sin(2.0*a/3.0 + M_PI/3);
}

static double fndd(double x, double y, double& dx, double& dy)
{
  double t1 = 2.0/3.0*atan2(x, y) + M_PI/3;
  double t2 = pow(x*x + y*y, 1.0/3.0);
  double t3 = x*x * ((y*y)/(x*x) + 1);
  dx = 2.0/3.0*x*sin(t1)/(t2*t2) + 2.0/3.0*y*t2*cos(t1)/t3;
  dy = 2.0/3.0*y*sin(t1)/(t2*t2) - 2.0/3.0*x*t2*cos(t1)/t3;
  return fn(x, y);
}


scalar bc_values(int marker, double x, double y)
{
  return fn(x, y);
}

scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}


int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("lshape.mesh");
  
  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
  space.set_bc_values(bc_values);
  space.set_uniform_order(2);
  
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form, SYM);
  
  ScalarView sview("Coarse solution", 0, 100, 798, 700);
  OrderView  oview("Polynomial orders", 800, 100, 798, 700);

  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Error Convergence for the L-shape Problem", "Degrees of Freedom", "Error [%]");
  graph.add_row("exact error", "k", "-", "o");
  graph.add_row("error estimate", "k", "--");

  Solution sln, rsln;
  UmfpackSolver umfpack;

  int it = 1;
  begin_time();
  while (1)
  {
    info("\n---- it=%d ------------------------------------------------------------------\n", it++);

    space.assign_dofs();
    
    LinSystem ls(&wf, &umfpack);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln);
    
    sview.show(&sln);
    oview.show(&space);

    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &rsln);

    ExactSolution exact(&mesh, fndd);
    double error = h1_error(&sln, &exact) * 100;
    info("\nExact solution error: %g%%", error);

    H1OrthoHP hp(1, &space);
    double estim = hp.calc_error(&sln, &rsln) * 100;
    info("Estimate of error: %g%%", estim);
    
    graph.add_values(0, space.get_num_dofs(), error);
    graph.add_values(1, space.get_num_dofs(), estim);
    graph.save("convergence.gp");
    
    if (error < 0.001) break;
    hp.adapt(0.3);
  }

  verbose("\nTotal run time: %g sec", end_time());
  View::wait();
  return 0;
}
