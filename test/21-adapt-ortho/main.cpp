#include "hermes2d.h"
//#include "solver_umfpack.h"

//#define LSHAPE


static double fn(double x, double y)
{
#ifndef LSHAPE
  return atan(60 * (sqrt(sqr(x-1.25) + sqr(y+0.25)) - M_PI/3));
#else
  double r = sqrt(x*x + y*y);
  double a = atan2(x, y);
  return pow(r, 2.0/3.0) * sin(2.0*a/3.0 + M_PI/3);
#endif
}

static double fndd(double x, double y, double& dx, double& dy)
{
#ifndef LSHAPE
  double t = sqrt(sqr(x-1.25) + sqr(y+0.25));
  double u = t * (3600 * sqr(t - M_PI/3) + 1);
  dx = 60 * (x-1.25) / u;
  dy = 60 * (y+0.25) / u;
  return fn(x, y);
#else
  double t1 = 2.0/3.0*atan2(x, y) + M_PI/3;
  double t2 = pow(x*x + y*y, 1.0/3.0);
  double t3 = x*x * ((y*y)/(x*x) + 1);
  dx = 2.0/3.0*x*sin(t1)/(t2*t2) + 2.0/3.0*y*t2*cos(t1)/t3;
  dy = 2.0/3.0*y*sin(t1)/(t2*t2) - 2.0/3.0*x*t2*cos(t1)/t3;
  return fn(x, y);
#endif
}


scalar bc_values(int marker, double x, double y)
{  
  return fn(x, y);
}

scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar bilinear_form_error(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv) + int_u_v(fu, fv, ru, rv);
}

double rhs(double x, double y)
{
  double t1 = sqrt(16*(x*x + y*y) - 40*x + 8*y + 26);
  double t2 = 3600*(x*x + y*y) - 9000*x + 1800*y;
  return -(240 * (t2 + 5849 - 400*M_PI*M_PI)) / (t1 * sqr(5851 + t2 - 600*t1*M_PI + 400*M_PI*M_PI));
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return -int_F_v(rhs, fv, rv);
}



int main(int argc, char* argv[])
{
  Mesh mesh, rmesh;
  #ifndef LSHAPE
  mesh.load("square.mesh");
  mesh.refine_all_elements();
  #else
  mesh.load("lshape.mesh");
  #endif
  
  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
  space.set_bc_values(bc_values);
  space.set_uniform_order(3);

  H1Space rspace(&rmesh, &shapeset);
  rspace.set_bc_values(bc_values);

  DiscreteProblem dp;
  dp.set_num_equations(1);
  dp.set_spaces(1, &space);
  dp.set_pss(1, &pss);
  dp.set_bilinear_form(0, 0, NULL, bilinear_form);
  #ifndef LSHAPE
  dp.set_linear_form(0, linear_form);
  #endif
  
  DiscreteProblem rp;
  rp.copy(&dp);
  rp.set_spaces(1, &rspace);

  ScalarView view("Solution", 0, 100, 798, 700);
  ScalarView ref("Refsol", 800, 100, 798, 700);
  OrderView  ord("Polynomial Orders", 800, 100, 798, 700);
  view.show_scale();
  view.show_mesh(false);

  GnuplotGraph graph;
  graph.set_log_y();
  graph.add_row("hp-ortho", "k", "-", "o");
  graph.add_row("Demkowicz", "k", "--");
  #ifndef LSHAPE
  graph.set_captions("Error Convergence for the Inner Layer Problem", "Degrees of Freedom", "Error [%]");
  double demkdof[] = { 50, 168, 415, 442, 666, 675, 1328, 1577, 1934, 2046, 2401, 2407, 2961 };
  double demkerr[] = { 77.7, 36.7, 18.9, 16.7, 12.3, 11.0, 5.76, 4.26, 1.78, 1.52, 1.185, 1.124, 0.6 };
  #else
  graph.set_captions("Error Convergence for the L-shape Problem", "Degrees of Freedom", "Error [%]");
  double demkdof[] = { 8, 22, 40, 107, 171, 234, 357, 419, 485, 549, 721, 837, 959 };
  double demkerr[] = { 22.5, 10.2, 6.35, 4.09, 2.57, 1.65, 0.755, 0.576, 0.483, 0.440, 0.155, 0.108, 0.07 };
  #endif
  graph.add_values(1, 13, demkdof, demkerr);

  int it = 1;
  bool done = false;
  begin_time();
  do
  {
    printf("\n\n---- it=%d ------------------------------------------------------------------\n\n", it++);

    space.assign_dofs();
    dp.create_matrix();
    dp.assemble_matrix_and_rhs();
    Solution sln;
    dp.solve_system(1, &sln);    

    view.show(&sln);
    ord.show(&space);
    //ord.wait_for_keypress();
    info("");

    rmesh.copy(&mesh);
    rmesh.refine_all_elements();
    rspace.copy_orders(&space, 1);
    rspace.assign_dofs();
    rp.create_matrix();
    rp.assemble_matrix_and_rhs();
    Solution rsln;
    rp.solve_system(1, &rsln);

    double err2 = 100 * h1_error(&sln, &rsln);
    info("\nNew Solution error: %g%%", err2);

    H1OrthoHP hp(1, &space);
    //double error = hp.calc_energy_error(&sln, &rsln, bilinear_form_error) * 100;
    double error = hp.calc_error(&sln, &rsln) * 100;
    if (error < 0.01) done = true;
    else hp.adapt(0.3, 1);

    info("\nEstimate of error: %g%%", error);
    graph.add_values(0, space.get_num_dofs(), error);
    #ifndef LSHAPE
    graph.save("conv_layer.gp");
    #else
    graph.save("conv_lshape.gp");
    #endif
    
    
/*    begin_time();
    H1OrthoHP hp(1, &space);
    double error = hp.calc_error(&sln, &rsln) * 100;
    info("\nError estimate: %g%%   (in %g sec)", error, end_time());
    #ifndef LSHAPE
    if (error < 0.005) break;
    #else
    if (error < 0.001) break;
    #endif
    
    begin_time();    
    hp.adapt(0.3);
    info("  (in %g sec)", end_time());*/
  }
  while (!done);

  verbose("\nTotal running time: %g sec", end_time());
  View::wait();
  return 0;
}

