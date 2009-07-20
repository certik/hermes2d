#include "hermes2d.h"
#include "solver_umfpack.h"

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

const int P_INIT = 1;             // initial polynomial degree in mesh
const double THRESHOLD = 0.3;     // the adaptivity algorithm goes on until THRESHOLD*total_error is processed
                                  // (see adapt_h1.cpp for explanation)
const int STRATEGY = 0;           // refinement strategy (0, 1, 2, 3 - see adapt_h1.cpp for explanation)
const bool H_ONLY = false;        // if H_ONLY == false then full hp-adaptivity takes place, otherwise
                                  // h-adaptivity is used. Use this parameter to check that indeed adaptive
                                  // hp-FEM converges much faster than adaptive h-FEM
const bool ISO_ONLY = false;      // when ISO_ONLY = true, only isotropic refinements are done,
                                  // otherwise also anisotropic refinements are allowed
const int MESH_REGULARITY = 1;    // specifies maximum allowed level of hanging nodes
                                  // -1 ... arbitrary level hangning nodes
                                  // 1, 2, 3,... means 1-irregular mesh, 2-irregular mesh, etc.
                                  // total regularization (0) is not supported in adaptivity
const double ERR_STOP = 0.001;    // adaptivity process stops when error wrt. exact solution in H1 norm
                                  // is less than this number
const int NDOF_STOP = 40000;      // adaptivity process stops when the number of degrees of freedom grows over
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
  graph.set_captions("Error Convergence for the L-Shape Domain Problem", "Degrees of Freedom", "Error [%]");
  graph.add_row("exact error", "k", "-", "o");
  graph.add_row("error estimate", "k", "--");

  GnuplotGraph graph_cpu;
  graph_cpu.set_captions("Error Convergence for the L-Shape Domain Problem", "CPU Time", "Error Estimate [%]");
  graph_cpu.add_row("exact error", "k", "-", "o");
  graph_cpu.add_row("error estimate", "k", "--");
  graph_cpu.set_log_y();

  Solution sln, rsln;
  UmfpackSolver umfpack;

  int it = 1;
  bool done = false;
  double cpu = 0.0;
  do
  {
    info("\n---- it=%d ------------------------------------------------------------------\n", it++);
    begin_time();

    // enumerating basis functions
    space.assign_dofs();

    LinSystem ls(&wf, &umfpack);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln);

    cpu += end_time();

    // view the solution
    sview.show(&sln);
    oview.show(&space);

    // solve the fine (reference) problem
    begin_time();
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &rsln);

    // calculate error wrt. exact solution
    ExactSolution exact(&mesh, fndd);
    double error = h1_error(&sln, &exact) * 100;
    info("\nExact solution error: %g%%", error);

    H1OrthoHP hp(1, &space);
    double err_est = hp.calc_error(&sln, &rsln) * 100;
    info("Estimate of error: %g%%", err_est);

    if (err_est < ERR_STOP || ls.get_num_dofs() >= NDOF_STOP) done = true;
    hp.adapt(THRESHOLD, STRATEGY, H_ONLY, ISO_ONLY, MESH_REGULARITY);
    cpu += end_time();

    // plotting convergence wrt. number of dofs
    graph.add_values(0, space.get_num_dofs(), error);
    graph.add_values(1, space.get_num_dofs(), err_est);
    graph.save("conv_dof.gp");

    // plotting convergence wrt. cpu time
    graph_cpu.add_values(0, cpu, error);
    graph_cpu.add_values(1, cpu, err_est);
    graph_cpu.save("conv_cpu.gp");
  }
  while (done == false);
  verbose("\nTotal running time: %g sec", cpu);

  View::wait();
  return 0;
}
