//  The purpose of this example is to show how to use Trilinos
//  for nonlinear PDE problem. It compares performance of NonlinSystem (Newton) with Umfpack
//  and performance of NOX solver (either using Newton's method or JFNK, with or without preconditioning)
//
//  PDE:  - \nabla (k \nabla u) = f
//  k = (1 + sqr(u_x) + sqr(u_y))^{-0.5}
//
//  Domain: square
//
//  BC: Dirichlet
//
//  Exact solution: (x - x*x) * (y - y*y)
//

#include "hermes2d.h"
#include "solver_umfpack.h"

const bool jfnk = false;     // true = jacobian-free method,
                            // false = Newton
const int precond = 2;      // preconditioning by jacobian (1) or approximation of jacobian (2)
                            // in case of jfnk,
                            // default ML proconditioner in case of Newton
const int ORDER = 3;

int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

///////////////////////////////////////////////////////////////////////////////////////////////
double exact(double x, double y, double &dx, double &dy)
{
	dx = (1- 2*x) * y * (1 - y);
	dy = (1- 2*y) * x * (1 - x);
	return  x * y * (1-x) * (1-y);
}

template<typename Real>
Real u(Real x, Real y)  {  return (x - x*x) * (y - y*y);  }
template<typename Real>
Real dudx(Real x, Real y)  {  return (1- 2*x) * y * (1 - y);  }
template<typename Real>
Real dudy(Real x, Real y)  {  return (1- 2*y) * x * (1 - x);  }

template<typename Real>
Real dudxx(Real x, Real y)  {  return -2.0 * (y-y*y);  }
template<typename Real>
Real dudyy(Real x, Real y)  {  return -2.0 * (x-x*x);  }
template<typename Real>
Real dudxy(Real x, Real y)  {  return (1- 2*y) * (1 - 2*x);  }

template<typename Real>
Real k(Real x, Real y)  {  return 1.0 / sqrt(1.0 + sqr(dudx(x,y)) + sqr(dudy(x,y)));  }
template<typename Real>
Real kx(Real x, Real y)  {  return -0.5 * pow(1.0 + sqr(dudx(x,y)) + sqr(dudy(x,y)), -1.5) *
                 (2.0 * dudx(x,y) * dudxx(x,y) + 2.0 * dudy(x,y) * dudxy(x,y));  }
template<typename Real>
Real ky(Real x, Real y)  {  return -0.5 * pow(1.0 + sqr(dudx(x,y)) + sqr(dudy(x,y)), -1.5) *
                 (2.0 * dudx(x,y) * dudxy(x,y) + 2.0 * dudy(x,y) * dudyy(x,y));  }

template<typename Real>
Real f(Real x, Real y)
{  return - kx(x,y) * dudx(x,y) - ky(x,y) * dudy(x,y) - k(x,y) * (dudxx(x,y) + dudyy(x,y)); }

///////////////////////////////////////////////////////////////////////////////////////////////
template<typename Real, typename Scalar>
Scalar jacobian_form_hermes(int n, double *wt, Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* u = ext->fn[0];
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( -0.5 * pow(1.0 + sqr(u->dx[i]) + sqr(u->dy[i]), -1.5) * (2.0 * u->dx[i] * vi->dx[i] + 2.0 * u->dx[i] * vi->dx[i])
                       * (u->dx[i] * vj->dx[i] + u->dy[i] * vj->dy[i]) +
                       (pow(1.0 + sqr(u->dx[i]) + sqr(u->dy[i]), -0.5))
                       * (vi->dx[i] * vj->dx[i] + vi->dy[i] * vj->dy[i]) );
  return result;
}

template<typename Real, typename Scalar>
Scalar residual_form_hermes(int n, double *wt, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* u = ext->fn[0];
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ((pow(1.0 + sqr(u->dx[i]) + sqr(u->dy[i]), -0.5)) * (u->dx[i] * vj->dx[i] + u->dy[i] * vj->dy[i])
                       - f(e->x[i], e->y[i]) * vj->val[i] );
  return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////
template<typename Real, typename Scalar>
Scalar jacobian_form_nox(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( -0.5 * pow(1.0 + sqr(u[0]->dx[i]) + sqr(u[0]->dy[i]), -1.5) *
                              (2.0 * u[0]->dx[i] * vi->dx[i] + 2.0 * u[0]->dx[i] * vi->dx[i])
                       * (u[0]->dx[i] * vj->dx[i] + u[0]->dy[i] * vj->dy[i]) +
                       (pow(1.0 + sqr(u[0]->dx[i]) + sqr(u[0]->dy[i]), -0.5))
                       * (vi->dx[i] * vj->dx[i] + vi->dy[i] * vj->dy[i]) );
  return result;
}

template<typename Real, typename Scalar>
Scalar precond_form_nox(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( vi->dx[i] * vj->dx[i] + vi->dy[i] * vj->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar residual_form_nox(int n, double *wt, Func<Real> *u[], Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ((pow(1.0 + sqr(u[0]->dx[i]) + sqr(u[0]->dy[i]), -0.5)) *
                        (u[0]->dx[i] * vj->dx[i] + u[0]->dy[i] * vj->dy[i])
                       - f(e->x[i], e->y[i]) * vj->val[i] );
  return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// main ////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  UmfpackSolver umfpack;
  Solution sln1, sln2;

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_uniform_order(ORDER);
  int ndofs = space.assign_dofs();

  info("******************** Using Linsystem, Solving by Umfpack ************************");

  begin_time();

  Solution prev;  prev.set_zero(&mesh);

  WeakForm wf1(1);
  wf1.add_biform(0, 0, callback(jacobian_form_hermes), UNSYM, ANY, 1, &prev);
  wf1.add_liform(0, callback(residual_form_hermes), ANY, 1, &prev);

  NonlinSystem sys(&wf1, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);
  sys.set_ic(&prev, &prev);

  bool done = false; int nit = 0;
  do
  {
    info("************* Newton: iteration %d ************", nit++);

    sys.assemble();
    sys.solve(1, &sln1);

    double res_l2_norm = sys.get_residuum_l2_norm();
    info("Residuum L2 norm: %g\n", res_l2_norm);
    if (res_l2_norm < 1e-6) done = true;

    prev.copy(&sln1);
  }
  while (!done);

  double umf_time = end_time();

  info("******************** Using FeProblem, Solving by NOX ************************");

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
  if (!jfnk || (jfnk && precond == 1)) wf2.add_jacform(0, 0, callback(jacobian_form_nox), SYM);
  if (jfnk && precond == 2) wf2.add_jacform(0, 0, callback(precond_form_nox), SYM);
  wf2.add_resform(0, callback(residual_form_nox));

  FeProblem fep(&wf2);
  fep.set_spaces(1, &space);
  fep.set_pss(1, &pss);

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
    vec = solver.get_solution();
    sln2.set_fe_solution(&space, &pss, vec);

    info("Number of nonlin iters: %d (norm of residual: %g)", solver.get_num_iters(), solver.get_residual());
    info("Total number of iters in linsolver: %d (achieved tolerance in the last step: %g)", solver.get_num_lin_iters(), solver.get_achieved_tol());
  }
  else
    error("Failed.");

  double nox_time = end_time();


  Solution ex;
  ex.set_exact(&mesh, &exact);
  info("\nSolution 1 (Hermes-Newton): h1 error: %g (time %g s)", 100 * h1_error(&sln1, &ex), umf_time);
  info("Solution 2 (Trilinos):      h1 error: %g (time %g + %g s)\n", 100 * h1_error(&sln2, &ex), proj_time, nox_time);

  ScalarView view1("Solution 1", 0, 0, 500, 400);
  view1.show(&sln1);

  ScalarView view2("Solution 2", 600, 0, 500, 400);
  view2.show(&sln2);




  View::wait("Waiting for all views to be closed.");
  return 0;
}
