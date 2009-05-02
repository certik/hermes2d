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
//  This is a very simple but nice example that allows you to compare h- and hp-adaptivity 
//  from the point of view of both CPU time requirements and discrete problem size, see
//  that the a-posteriori error estimate Hermes uses is almost as good as the exact error, etc. 
//  The following problem parameters can be changed easily:
//

int P_INIT = 1;           // initial polynomial degree in mesh
double THRESHOLD = 0.3;   // the adaptivity algorithm goes on until THRESHOLD*total_error is processed
                          // (see adapt_h1.cpp for explanation)
int STRATEGY = 0;         // refinement strategy (0, 1, 2, 3 - see adapt_h1.cpp for explanation)
int H_ONLY = 0;           // if H_ONLY == 0 then full hp-adaptivity takes place, otherwise
                          // h-adaptivity is used. Use this parameter to check that indeed adaptive 
                          // hp-FEM converges much faster than adaptive h-FEM
double ERR_STOP = 0.001;  // adaptivity process stops when error wrt. exact solution in H1 norm 
                          // is less than this number 
int NDOF_STOP = 40000;    // adaptivity process stops when the number of degrees of freedom grows over 
                          // this limit. This is mainly to prevent h-adaptivity to go on forever.  

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
  if(P_INIT == 1) mesh.refine_all_elements();  // this is because there are no degrees of freedom 
                                               // on the coarse mesh lshape.mesh if P_INIT == 1
  
  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
  
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
    
    if (error < ERR_STOP || ls.get_num_dofs() >= NDOF_STOP) break;
    hp.adapt(THRESHOLD, STRATEGY, H_ONLY);
  }

  verbose("\nTotal run time: %g sec", end_time());
  View::wait();
  return 0;
}
