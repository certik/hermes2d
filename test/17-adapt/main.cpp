#include "hermes2d.h"


static double fn(double x, double y)
{
  return atan(60 * (sqrt(sqr(x-1.25) + sqr(y+0.25)) - M_PI/3));
}

static double fndd(double x, double y, double& dx, double& dy)
{
  double t = sqrt(sqr(x-1.25) + sqr(y+0.25));
  double u = t * (3600 * sqr(t - M_PI/3) + 1);
  dx = 60 * (x-1.25) / u;
  dy = 60 * (y+0.25) / u;
  return fn(x, y);
}


double bc_values(int marker, double x, double y)
{  
  return fn(x, y);
}

double bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return -int_grad_u_grad_v(fu, fv, ru, rv);
}

static double rhs(double x, double y)
{
  double t1 = sqrt(16*(x*x + y*y) - 40*x + 8*y + 26);
  double t2 = 3600*(x*x + y*y) - 9000*x + 1800*y;
  return -(240 * (t2 + 5849 - 400*M_PI*M_PI)) / (t1 * sqr(5851 + t2 - 600*t1*M_PI + 400*M_PI*M_PI));
}

double linear_form(RealFunction* fv, RefMap* rv)
{
  return int_F_v(rhs, fv, rv);
}


int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);

  Mesh mesh;
  mesh.load("square.mesh");
  //mesh.load("curved_quad2.mesh");
  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
  space.set_bc_values(bc_values);
  space.set_uniform_order(1);

  DiscreteProblem dp;
  dp.set_num_equations(1);
  dp.set_spaces(1, &space);
  dp.set_pss(1, &pss);
  dp.set_bilinear_form(0, 0, bilinear_form);
  dp.set_linear_form(0, linear_form);

  ScalarView view("Solution", 0, 100, 798, 700);
  ScalarView ref("Refsol", 800, 100, 798, 798);
  OrderView  ord("Polynomial Orders", 800, 100, 798, 700);
  view.show_scale();
  view.show_mesh(false);

  MatlabGraph graph;
  graph.set_captions("Error Convergence for the Inner Layer Problem", "Degrees of Freedom", "Error [%]");
  graph.set_log_y();
  graph.add_row("Our Algorithm", "k", "-", "o");
  graph.add_row("Demkowicz et al.", "k", "--");
  graph.add_row("Estimated error", "k", ":");
  double demkdof[] = { 50, 168, 415, 442, 666, 675, 1328, 1577, 1934, 2046, 2401, 2407, 2961 };
  double demkerr[] = { 77.7, 36.7, 18.9, 16.7, 12.3, 11.0, 5.76, 4.26, 1.78, 1.52, 1.185, 1.124, 0.6 };
  graph.add_values(1, 13, demkdof, demkerr);

  int it = 1;
  bool cont;
  begin_time();
  do
  {
    printf("\n\n---- it=%d ------------------------------------------------------------------\n\n", it++);

    space.assign_dofs();
    dp.create_matrix();
    dp.assemble_matrix_and_rhs();

    Solution sln;
    dp.solve_system(1, &sln);    

    double error = 100 * h1_error_norm_exact(&sln, fndd);
    info("\nExact solution error: %g%%", error);
    graph.add_values(0, space.get_num_dofs(), error);
    
    view.show(&sln);
    ord.show(&space);

    double2 energy2[1];
    cont = adapt_energy(&dp, 0.1, energy2);
    graph.add_values(2, space.get_num_dofs(), sqrt(energy2[0][0] / energy2[0][1]) * 100);
    
    graph.save("conv.m");
  }
  while (cont);

  verbose("\nTotal running time: %g sec", end_time());
  hermes2d_finalize();
  return 0;
}
