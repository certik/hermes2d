//  The purpose of this example is to show how to use Trilinos
//  for time-dependent PDE problem.
//  NOX solver is used, either using Newton's method or JFNK and
//  with or without preconditioning
//
//  PDE: heat-transfer
//
//  Domain: square
//
//  BC: cooling Dirichlet at the bottom
//

#include "hermes2d.h"
#include "solver_umfpack.h"

const double ALPHA = 10.0;
const double LAMBDA = 1e5;
const double HEATCAP = 1e6;
const double RHO = 3000.0;
const double TEMP_EXT = 20.0;
const double TEMP_INIT = 10.0;

const double tau = 50.0;

const bool jfnk = true;
const bool precond = true;

const int ORDER = 1;

//////////////////////////////////////////////////////////////////////////////////////////////
int bc_types(int marker)
{
  if (marker == 1) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

scalar bc_values(int marker, double x, double y)
{
  return TEMP_INIT;
}

//////////////////////////////////////////////////////////////////////////////////////////////
template<typename Real, typename Scalar>
Scalar jacobian(int n, double *wt, Func<Scalar> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (HEATCAP * RHO * vi->val[i] * vj->val[i] / tau
                     + LAMBDA * (vi->dx[i] * vj->dx[i] + vi->dy[i] * vj->dy[i]));
  return result;
}

template<typename Real, typename Scalar>
Scalar jacobian_surf(int n, double *wt, Func<Scalar> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] *(LAMBDA * ALPHA * vi->val[i] * vj->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar residual(int n, double *wt, Func<Scalar> *u[], Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (HEATCAP * RHO * (u[0]->val[i] - ext->fn[0]->val[i]) * vj->val[i] / tau
                     + LAMBDA * (u[0]->dx[i] * vj->dx[i] + u[0]->dy[i] * vj->dy[i]));
  return result;
}

template<typename Real, typename Scalar>
Scalar residual_surf(int n, double *wt, Func<Scalar> *u[], Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (LAMBDA * ALPHA * (u[0]->val[i] - TEMP_EXT) * vj->val[i]);
  return result;
}

////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  // load and refine mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  // set up shapeset
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // set up spaces
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(ORDER);
  int ndofs = space.assign_dofs();

  // set initial condition
  Solution tprev, titer, tsln;
  tprev.set_const(&mesh, 20.0);
  titer.set_const(&mesh, 20.0);

  WeakForm wf(1, jfnk ? true : false);
  wf.add_jacform(0, 0, callback(jacobian));
  wf.add_jacform_surf(0, 0, callback(jacobian_surf));
  wf.add_resform(0, callback(residual), ANY, 1, &tprev);
  wf.add_resform_surf(0, callback(residual_surf));

  // Finite element problem
  FeProblem fep(&wf);
  fep.set_spaces(1, &space);
  fep.set_pss(1, &pss);

  // Obtain solution vector for initial guess
  begin_time();
  info("Projecting initial solution");
  Projection proj(1, &titer, &space, &pss);
  UmfpackSolver umfpack;
  proj.set_solver(&umfpack);
  double* vec = proj.project();
  double proj_time = end_time();

  // Solver + preconditioner
  NoxSolver solver(&fep);
  MlPrecond pc("sa");
  if (precond)
  {
    if (jfnk) solver.set_precond(&pc);
    else solver.set_precond("ML");
  }

  // visualisation
  ScalarView Tview("Temperature", 0, 0, 450, 600);
  Tview.set_min_max_range(10,20);

  double total_time = 0.0;
  begin_time();
  for (int it = 1; total_time <= 2000.0; it++)
  {
    info("\n*** Time iteration %d, t = %g s ***", it, total_time += tau);

    solver.set_init_sln(vec);
    bool solved = solver.solve();
    if (solved)
    {
      vec = solver.get_solution();
      tsln.set_fe_solution(&space, &pss, vec);
      Tview.show(&tsln);
      tprev = tsln;
    }
    else
      error("Failed.");

    info("Number of nonlin iters: %d (norm of residual: %g)", solver.get_num_iters(), solver.get_residual());
    info("Total number of iters in linsolver: %d (achieved tolerance in the last step: %g)", solver.get_num_lin_iters(), solver.get_achieved_tol());
  }

  info("Total running time: %g", end_time());

  View::wait("Waiting for all views to be closed.");
  return 0;
}
