//  The purpose of this example is to show how to use Trilinos for nonlinear time-dependent coupled system of PDEs
//  Solved by NOX solver, JFNK (derivation of Jacobian is not needed), with or without preconditioning
//
//  PDE: Flame propagation
//
//  Domain: pipe with rods
//

#include "hermes2d.h"
#include "solver_umfpack.h"

const bool jfnk = false;     // true = jacobian-free method,
                            // false = Newton
const int precond = 2;      // preconditioning by jacobian (1) (less GMRES iterations, more time to create precond)
                            // or by approximation of jacobian (2) (less time for precond creation, more GMRES iters)
                            // in case of jfnk,
                            // default Ifpack proconditioner in case of Newton

const bool trilinos_output = false;  // display more details about nonlinear and linear solvers


// problem constants - all according to the paper SchmichVexler2008
const double Le    = 1.0;
const double alpha = 0.8;
const double beta  = 10.0;
const double kappa = 0.1;
const double x1    = 9.0;

const double tau = 0.05;
const int init_order = 2;

//// BCs, omega ////////////////////////////////////////////////////////////////////////////////////

int bc_types(int marker)
  { return (marker == 1) ? BC_ESSENTIAL : BC_NATURAL; }

scalar temp_bc_values(int marker, double x, double y)
  { return (marker == 1) ? 1.0 : 0; }


scalar temp_ic(double x, double y, scalar& dx, scalar& dy)
  { return (x <= x1) ? 1.0 : exp(x1 - x); }

scalar conc_ic(double x, double y, scalar& dx, scalar& dy)
  { return (x <= x1) ? 0.0 : 1.0 - exp(Le*(x1 - x)); }

////////////////////////////////////////////////////////////////////////////////////////////////////////
void omega_fn(int n, scalar* a, scalar* dadx, scalar* dady, scalar* b, scalar* dbdx, scalar* dbdy,
                      scalar* out, scalar* outdx, scalar* outdy)
{
  for (int i = 0; i < n; i++)
  {
    double t1 = a[i] - 1.0;
    double t2 = t1 * beta;
    double t3 = 1.0 + t1 * alpha;
    double t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    double t5 = (beta / (t3 * t3)) * b[i];
    out[i] = t4 * b[i];
    outdx[i] = t4 * (dbdx[i] + dadx[i] * t5);
    outdy[i] = t4 * (dbdy[i] + dady[i] * t5);
  }
}

void omega_dt_fn(int n, scalar* a, scalar* dadx, scalar* dady, scalar* b, scalar* dbdx, scalar* dbdy,
                        scalar* out, scalar* outdx, scalar* outdy)
{
  for (int i = 0; i < n; i++)
  {
    double t1 = a[i] - 1.0;
    double t2 = t1 * beta;
    double t3 = 1.0 + t1 * alpha;
    double t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    double t5 = (beta / (t3 * t3));
    out[i] = t4 * t5 * b[i];
    outdx[i] = 0.0;
    outdy[i] = 0.0; // not important
  }
}

void omega_dc_fn(int n, scalar* a, scalar* dadx, scalar* dady, scalar* b, scalar* dbdx, scalar* dbdy,
                        scalar* out, scalar* outdx, scalar* outdy)
{
  for (int i = 0; i < n; i++)
  {
    double t1 = a[i] - 1.0;
    double t2 = t1 * beta;
    double t3 = 1.0 + t1 * alpha;
    double t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    out[i] = t4;
    outdx[i] = 0.0;
    outdy[i] = 0.0; // not important
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename Real>
Real omega(Real t, Real c)
{
  Real t2 = beta * (t - 1.0);
  Real t3 = 1.0 + alpha * (t - 1.0);
  return c * sqr(beta) / (2.0*Le) * exp(t2 / t3);
}

/////////   Residual   ///////////////////////////////////////////////////////////////////////////////////////////
template<typename Real, typename Scalar>
Scalar residual_0(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* tprev1 = ext->fn[0];
  Func<Real>* tprev2 = ext->fn[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( (3.0 * u[0]->val[i] - 4.0 * tprev1->val[i] + tprev2->val[i]) * vi->val[i] / (2.0 * tau) +
                        (u[0]->dx[i] * vi->dx[i] + u[0]->dy[i] * vi->dy[i]) -
                        omega(u[0]->val[i], u[1]->val[i]) * vi->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar residual_0_surf(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (kappa * u[0]->val[i] * vi->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar residual_1(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* cprev1 = ext->fn[0];
  Func<Real>* cprev2 = ext->fn[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( (3.0 * u[1]->val[i] - 4.0 * cprev1->val[i] + cprev2->val[i]) * vi->val[i] / (2.0 * tau) +
                        (u[1]->dx[i] * vi->dx[i] + u[1]->dy[i] * vi->dy[i]) / Le +
                        omega(u[0]->val[i], u[1]->val[i]) * vi->val[i] );
  return result;
}

//////////   Preconditioning   /////////////////////////////////////////////////////////////////////////////////////////////
template<typename Real, typename Scalar>
Scalar precond_0_0(int n, double *wt, Func<Scalar>* u[], Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / tau
                      +  vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar precond_1_1(int n, double *wt,  Func<Scalar>* u[], Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / tau
                      +  (vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]) / Le );
  return result;
}

//////   Jacobian   /////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Real, typename Scalar>
Scalar jacobian_0_0(int n, double *wt, Func<Scalar>* u[], Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / tau
                      +  vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]
                      - (u[1]->val[i] * sqr(beta)/(2*Le) *
                         exp(beta * (u[0]->val[i] - 1)/(1.0 + alpha*(u[0]->val[i] - 1))) *
                         beta / (sqr(1.0 + alpha*(u[0]->val[i] - 1))))
                           * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar jacobian_0_0_surf(int n, double *wt, Func<Scalar>* u[], Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (kappa * vj->val[i] * vi->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar jacobian_0_1(int n, double *wt, Func<Scalar>* u[], Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* dodc = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (- (sqr(beta)/(2*Le) * exp(beta * (u[0]->val[i]-1)/(1.0 + alpha*(u[0]->val[i]-1)))) * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar jacobian_1_0(int n, double *wt, Func<Scalar>* u[], Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* dodt = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( (u[1]->val[i] * sqr(beta)/(2*Le) *
                         exp(beta * (u[0]->val[i] - 1)/(1.0 + alpha*(u[0]->val[i] - 1))) *
                         beta / (sqr(1.0 + alpha*(u[0]->val[i] - 1)))) * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar jacobian_1_1(int n, double *wt,  Func<Scalar>* u[], Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / tau
                      +  (vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]) / Le
                      + (sqr(beta)/(2*Le) * exp(beta * (u[0]->val[i]-1)/(1.0 + alpha*(u[0]->val[i]-1))))
                        * vj->val[i] * vi->val[i] );
  return result;
}

//// main //////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  H1Shapeset shapeset;
  PrecalcShapeset tpss(&shapeset);
  PrecalcShapeset cpss(&shapeset);

  // initialize meshes
  H1Space tspace(&mesh, &shapeset);
  H1Space cspace(&mesh, &shapeset);
  tspace.set_bc_types(bc_types);
  tspace.set_bc_values(temp_bc_values);
  cspace.set_bc_types(bc_types);
  tspace.set_uniform_order(init_order);
  cspace.set_uniform_order(init_order);
  int ndofs = 0;
  ndofs += tspace.assign_dofs(ndofs);
  ndofs += cspace.assign_dofs(ndofs);
  info("Number of DOFs: %d", ndofs);

  // initial conditions
  Solution tprev1, cprev1, tprev2, cprev2, titer, citer, tsln, csln;
  tprev1.set_exact(&mesh, temp_ic);  cprev1.set_exact(&mesh, conc_ic);
  tprev2.set_exact(&mesh, temp_ic);  cprev2.set_exact(&mesh, conc_ic);
  titer.set_exact(&mesh, temp_ic);   citer.set_exact(&mesh, conc_ic);
  DXDYFilter omega_dt(omega_dt_fn, &tprev1, &cprev1);
  DXDYFilter omega_dc(omega_dc_fn, &tprev1, &cprev1);

  begin_time();
  info("Projecting initial solution");
  Projection proj(2, &titer, &citer, &tspace, &cspace, &tpss, &cpss);
  UmfpackSolver umfpack;
  proj.set_solver(&umfpack);
  double* vec = proj.project();
  double proj_time = end_time();

  //   visualization
  ScalarView rview("Reaction rate", 0, 0, 1600, 460);
  rview.set_min_max_range(0.0,2.0);

  WeakForm wf(2, jfnk ? true : false);
  if (!jfnk || (jfnk && precond == 1))
  {
    wf.add_jacform(0, 0, callback(jacobian_0_0));
    wf.add_jacform_surf(0, 0, callback(jacobian_0_0_surf));
    wf.add_jacform(1, 1, callback(jacobian_1_1));
    wf.add_jacform(0, 1, callback(jacobian_0_1));
    wf.add_jacform(1, 0, callback(jacobian_1_0));
  }
  else if (precond == 2)
  {
    wf.add_jacform(0, 0, callback(precond_0_0));
    wf.add_jacform(1, 1, callback(precond_1_1));
  }

  wf.add_resform(0, callback(residual_0), ANY, 2, &tprev1, &tprev2);
  wf.add_resform_surf(0, callback(residual_0_surf), 3);
  wf.add_resform(1, callback(residual_1), ANY, 2, &cprev1, &cprev2);

  FeProblem fep(&wf);
  fep.set_spaces(2, &tspace, &cspace);
  fep.set_pss(2, &tpss, &cpss);

  NoxSolver solver(&fep);
  MlPrecond pc("sa");
  if (precond)
  {
    if (jfnk) solver.set_precond(&pc);
    else solver.set_precond("Ifpack");
  }
  if (trilinos_output)
    solver.set_output_flags(NOX::Utils::Error | NOX::Utils::OuterIteration |
                            NOX::Utils::OuterIterationStatusTest |
                            NOX::Utils::LinearSolverDetails);

  double total_time = 0.0;
  for (int it = 1; total_time <= 60.0; it++)
  {
    info("\n*** Time iteration %d, t = %g s ***", it, total_time + tau);

    begin_time();
    solver.set_init_sln(vec);
    bool solved = solver.solve();
    if (solved)
    {
      vec = solver.get_solution();
      tsln.set_fe_solution(&tspace, &tpss, vec);
      csln.set_fe_solution(&cspace, &cpss, vec);

      info("Number of nonlin iters: %d (norm of residual: %g)",
          solver.get_num_iters(), solver.get_residual());
      info("Total number of iters in linsolver: %d (achieved tolerance in the last step: %g)",
          solver.get_num_lin_iters(), solver.get_achieved_tol());

      // visualization
      DXDYFilter omega_view(omega_fn, &tsln, &csln);
      rview.show(&omega_view);

      total_time += tau;

      tprev2.copy(&tprev1);
      cprev2.copy(&cprev1);
      tprev1 = tsln;
      cprev1 = csln;
    }
    else
      error("Failed.");

    info("Total running time for time level %d: %g s.", it, end_time());
  }

  View::wait("Waiting for all views to be closed.");
  return 0;
}

