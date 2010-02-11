//  The purpose of this example is to show how to use Trilinos while adapting mesh
//  Solved by NOX solver, either using Newton's method or JFNK, with or without preconditioning
//
//  PDE: -Laplace u = f
//
//  Domain: square
//
//  BC: Dirichlet
//
//  Known exact solution, see functions fn() and fndd()
//

#include "hermes2d.h"
#include "solver_umfpack.h"

const bool jfnk = true;
const bool precond = true;

const int INIT_ORDER = 2;
const double ERR_STOP = 1.0;

// problem constants
double SLOPE = 200;       // slope of the step inside the domain

// exact solution
static double fn(double x, double y)
{
  return atan(SLOPE * (sqrt(sqr(x-1.25) + sqr(y+0.25)) - M_PI/3));
}

static double fndd(double x, double y, double& dx, double& dy)
{
  double t = sqrt(sqr(x-1.25) + sqr(y+0.25));
  double u = t * (sqr(SLOPE) * sqr(t - M_PI/3) + 1);
  dx = SLOPE * (x-1.25) / u;
  dy = SLOPE * (y+0.25) / u;
  return fn(x, y);
}

int bc_types(int marker)
{
  return BC_ESSENTIAL;
}
// boundary conditions
scalar bc_values(int marker, double x, double y)
{
  return fn(x, y);
}


template<typename Real>
Real rhs(Real x, Real y)
{
  Real t2 = sqr(y + 0.25) + sqr(x - 1.25);
  Real t = sqrt(t2);
  Real u = (sqr(M_PI - 3.0*t)*sqr(SLOPE) + 9.0);
  return 27.0/2.0 * sqr(2.0*y + 0.5) * (M_PI - 3.0*t) * pow(SLOPE,3.0) / (sqr(u) * t2) +
         27.0/2.0 * sqr(2.0*x - 2.5) * (M_PI - 3.0*t) * pow(SLOPE,3.0) / (sqr(u) * t2) -
          9.0/4.0 * sqr(2.0*y + 0.5) * SLOPE / (u * pow(t,3.0)) -
          9.0/4.0 * sqr(2.0*x - 2.5) * SLOPE / (u * pow(t,3.0)) +
          18.0 * SLOPE / (u * t);
}

template<typename Real, typename Scalar>
Scalar precond_form(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, vi, vj);
}

template<typename Real, typename Scalar>
Scalar residual_form(int n, double *wt, Func<Real> *u[], Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u[0], vj) + int_F_v<Real, Scalar>(n, wt, rhs, vj, e);
}


int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);
  mesh.refine_all_elements();

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create finite element space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(INIT_ORDER);

  // initialize the weak formulation
  WeakForm wf(1, jfnk ? true : false);
  if (precond) wf.add_jacform(0, 0, callback(precond_form), SYM);
  wf.add_resform(0, callback(residual_form));

  // visualize solution and mesh
  ScalarView sview("Coarse solution", 0, 100, 798, 700);
  OrderView  oview("Polynomial orders", 800, 100, 798, 700);

  // convergence graph wrt. the number of degrees of freedom
  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Error Convergence for the Inner Layer Problem", "Degrees of Freedom", "Error [%]");
  graph.add_row("exact error", "k", "-", "o");
  graph.add_row("error estimate", "k", "--");

  // adaptivity loop
  int it = 0, ndofs;
  bool done = false;
  double error;
  Solution sln_coarse, sln_fine;
  double *s;
  do
  {
    info("\n---- Adaptivity step %d ---------------------------------------------\n", it++);

    ndofs = space.assign_dofs();
    FeProblem fep(&wf);
    fep.set_spaces(1, &space);
    fep.set_pss(1, &pss);
    NoxSolver solver(&fep);
    MlPrecond pc("sa");
    if (precond)
    {
      if (jfnk) solver.set_precond(&pc);
      else solver.set_precond("ML");
    }
    bool solved = solver.solve();
    if (solved)
    {
      s = solver.get_solution();
      sln_coarse.set_fe_solution(&space, &pss, s);

      info("\nCoarse Solution info:");
      info("Number of nonlin iters: %d (norm of residual: %g)", solver.get_num_iters(), solver.get_residual());
      info("Total number of iters in linsolver: %d (achieved tolerance in the last step: %g)", solver.get_num_lin_iters(), solver.get_achieved_tol());

      // view the solution and mesh
      sview.show(&sln_coarse);
      oview.show(&space);
    }

    Mesh rmesh; rmesh.copy(&mesh); rmesh.refine_all_elements();
    H1Space rspace(&rmesh, &shapeset);
    rspace.set_bc_types(bc_types);
    rspace.set_bc_values(bc_values);
    rspace.copy_orders(&space, 1);
    rspace.assign_dofs();

    FeProblem ref_fep(&wf);
    ref_fep.set_spaces(1, &rspace);
    ref_fep.set_pss(1, &pss);
    NoxSolver ref_solver(&ref_fep);
    if (precond)
    {
      if (jfnk) ref_solver.set_precond(&pc);
      else ref_solver.set_precond("ML");
    }
    solved = ref_solver.solve();
    if (solved)
    {
      s = ref_solver.get_solution();
      sln_fine.set_fe_solution(&rspace, &pss, s);

      info("\nReference Solution info:");
      info("Number of nonlin iters: %d (norm of residual: %g)",
            ref_solver.get_num_iters(), ref_solver.get_residual());
      info("Total number of iters in linsolver: %d (achieved tolerance in the last step: %g)",
            ref_solver.get_num_lin_iters(), ref_solver.get_achieved_tol());
    }
    else
      error("Failed.");

    // calculate error estimate wrt. fine mesh solution
    H1OrthoHP hp(1, &space);
    double err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;
    info("\nEstimate of error: %g%%  (ndofs %d)", err_est, ndofs);
    ExactSolution exact(&mesh, fndd);
    error = h1_error(&sln_coarse, &exact) * 100;
    info("Exact solution error: %g%%", error);

    // add entry to DOF convergence graph
    graph.add_values(0, space.get_num_dofs(), error);
    graph.add_values(1, space.get_num_dofs(), err_est);
    graph.save("conv_dof.gp");

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else  hp.adapt(0.3, 1, 0);
  }
  while (done == false);

  View::wait("Waiting for all views to be closed.");
  return 0;
}
