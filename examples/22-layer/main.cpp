#include "hermes2d.h"
#include "solver_umfpack.h"

// 
//  PDE: -Laplace u = 0
//
//  Known exact solution, see functions fn() and fndd() 
//
//  Domain: L-shape domain, see the file lshape.mesh
//
//  BC:  Dirichlet, given by exact solution
//
//  In addition to the L-shape domain problem, this is another example that allows you to compare 
//  h- and hp-adaptivity from the point of view of both CPU time requirements and discrete problem 
//  size, look at the quality of the a-posteriori error estimator used by Hermes (exact error is also)
//  provided, etc. The following problem parameters can be changed easily:
//

int P_INIT = 1;           // initial polynomial degree in mesh
double THRESHOLD = 0.3;   // the adaptivity algorithm goes on until THRESHOLD*total_error is processed
                          // (see adapt_h1.cpp for explanation)
int STRATEGY = 0;         // refinement strategy (0, 1, 2, 3 - see adapt_h1.cpp for explanation)
int H_ONLY = 0;           // if H_ONLY == 0 then full hp-adaptivity takes place, otherwise
                          // h-adaptivity is used. Use this parameter to check that indeed adaptive 
                          // hp-FEM converges much faster than adaptive h-FEM
double ERR_STOP = 0.01;  // adaptivity process stops when error wrt. exact solution in H1 norm 
                          // is less than this number 
int NDOF_STOP = 40000;    // adaptivity process stops when the number of degrees of freedom grows over 
                          // this limit. This is mainly to prevent h-adaptivity to go on forever.  

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


scalar bc_values(int marker, double x, double y)
{  
  return fn(x, y);
}

scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
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
  Mesh mesh;
  mesh.load("square_quad.mesh");
  if(P_INIT == 1) mesh.refine_all_elements();  // this is because there are no degrees of freedom 
                                               // on the coarse mesh lshape.mesh if P_INIT == 1
 
  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
  
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form, SYM);
  wf.add_liform(0, linear_form);

  ScalarView sview("Coarse solution", 0, 100, 798, 700);
  OrderView  oview("Polynomial orders", 800, 100, 798, 700);

  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Error Convergence for the Inner Layer Problem", "Degrees of Freedom", "Error [%]");
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
    
    if (error < ERR_STOP || ls.get_num_dofs() >= NDOF_STOP) break;
    hp.adapt(THRESHOLD, STRATEGY, H_ONLY);
  }

  verbose("\nTotal run time: %g sec", end_time());
  View::wait();
  return 0;
}
