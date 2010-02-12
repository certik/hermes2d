//  The purpose of this example is to show how to use Trilinos
//  for linear PDE problem. It compares performance of LinSystem with Umfpack
//  and performance of NOX solver (either using Newton's method or JFNK,
//  with or without preconditioning)
//
//  PDE: Poisson equation
//
//  Domain: square
//
//  BC: Dirichlet
//
//  Exact solution: sqr(x) + sqr(y)
//

#include "hermes2d.h"
#include "solver_umfpack.h"

const bool jfnk = true;     // true = jacobian-free method,
                            // false = Newton
const bool precond = true;  // preconditioning by jacobian in case of jfnk,
                            // default ML proconditioner in case of Newton
const int ORDER = 1;

int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

double exact(double x, double y, double &dx, double &dy)
{
	dx = 2*x;
	dy = 2*y;
	return x*x +y*y;
}

double bc_values(int marker, double x, double y)
{
  return x*x + y*y;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ((-4.0) * v->val[i]);
  return result;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Real, typename Scalar>
Scalar jacobian_form(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( vi->dx[i] * vj->dx[i] + vi->dy[i] * vj->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar residual_form(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( u[0]->dx[i] * vi->dx[i] + u[0]->dy[i] * vi->dy[i] + 4.0 * vi->val[i] );
  return result;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(ORDER);
  int ndofs = space.assign_dofs();

  info("\n******************** Using Linsystem, Solving by Umfpack ************************");

  begin_time();

  UmfpackSolver umfpack;
  Solution sln1;

  WeakForm wf1(1);
  wf1.add_biform(0, 0, callback(bilinear_form));
  wf1.add_liform(0, callback(linear_form));

  LinSystem sys(&wf1, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);

  sys.assemble();
  sys.solve(1, &sln1);

  double umf_time = end_time();

  info("\n******************** Using FeProblem, Solving by NOX Solver ************************");

  begin_time();
  info("Projecting initial solution");
  Solution init;  init.set_zero(&mesh);
  Projection proj(1, &init, &space, &pss);
  proj.set_solver(&umfpack);
  double* vec = proj.project();
  double proj_time = end_time();

  begin_time();
  info("Number of DOFs: %d", ndofs);
  WeakForm wf2(1, jfnk ? true : false);
  wf2.add_jacform(0, 0, callback(jacobian_form), SYM);
  wf2.add_resform(0, callback(residual_form));

  FeProblem fep(&wf2);
  fep.set_spaces(1, &space);
  fep.set_pss(1, &pss);

  Solution sln2;
  NoxSolver solver(&fep);
  solver.set_init_sln(vec);
  MlPrecond pc("sa");
  if (precond)
  {
    if (jfnk) solver.set_precond(&pc);
    else solver.set_precond("ML");
  }
  bool solved = solver.solve();
  if (solved)
  {
    double *s = solver.get_solution();
    sln2.set_fe_solution(&space, &pss, s);
    info("Number of nonlin iters: %d (norm of residual: %g)", solver.get_num_iters(), solver.get_residual());
    info("Total number of iters in linsolver: %d (achieved tolerance in the last step: %g)", solver.get_num_lin_iters(), solver.get_achieved_tol());

  }
  else
    error("Failed");

  double nox_time = end_time();

  Solution ex;
  ex.set_exact(&mesh, &exact);
  info("\nSolution 1 (Umfpack):  h1 error: %g (time %g s)", 100 * h1_error(&sln1, &ex), umf_time);
  info("Solution 2 (Trilinos): h1 error: %g (time %g + %g s)\n", 100 * h1_error(&sln2, &ex), proj_time, nox_time);

  ScalarView view1("Solution 1", 0, 0, 500, 400);
  view1.set_min_max_range(0, 2);
  view1.show(&sln1);

  ScalarView view2("Solution 2", 600, 0, 500, 400);
  view2.set_min_max_range(0, 2);
  view2.show(&sln2);


  View::wait("Waiting for all views to be closed.");
  return 0;
}
