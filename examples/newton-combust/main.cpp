#include "hermes2d.h"
#include "solver_umfpack.h"

//  This is a flame propagation problem
//
//  PDEs:
//
//  dT/dt - laplace T = omega(T,Y)
//  dY/dt - 1/Le * laplace Y = - omega(T,Y)
//
//  Domain: rectangle with cooled rods
//
//  BC:  T = 1, Y = 0 on the inlet
//       dT/dn = - kappa T on cooled rods
//       dT/dn = 0, dY/dn = 0 elsewhere
//
//  Time-stepping: second order BDF formula

// problem constants - all according to the paper SchmichVexler2008
const double Le    = 1.0;
const double alpha = 0.8;
const double beta  = 10.0;
const double kappa = 0.1;
const double x1    = 9.0;

const double tau = 0.05;
const int init_order = 3;

//// BCs, omega ////////////////////////////////////////////////////////////////////////////////////

int bc_types(int marker)
  { return (marker == 1) ? BC_ESSENTIAL : BC_NATURAL; }

scalar temp_bc_values(int marker, double x, double y)
  { return (marker == 1) ? 1.0 : 0; }


scalar temp_ic(double x, double y, scalar& dx, scalar& dy)
  { return (x <= x1) ? 1.0 : exp(x1 - x); }

scalar conc_ic(double x, double y, scalar& dx, scalar& dy)
  { return (x <= x1) ? 0.0 : 1.0 - exp(Le*(x1 - x)); }


void omega_fn(int n, scalar* a, scalar* dadx, scalar* dady, scalar* b, scalar* dbdx, scalar* dbdy,
                      scalar* out, scalar* outdx, scalar* outdy)
{
  for (int i = 0; i < n; i++)
  {
    scalar t1 = a[i] - 1.0;
    scalar t2 = t1 * beta;
    scalar t3 = 1.0 + t1 * alpha;
    scalar t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    scalar t5 = (beta / (t3 * t3)) * b[i];
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
    scalar t1 = a[i] - 1.0;
    scalar t2 = t1 * beta;
    scalar t3 = 1.0 + t1 * alpha;
    scalar t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    scalar t5 = (beta / (t3 * t3));
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
    scalar t1 = a[i] - 1.0;
    scalar t2 = t1 * beta;
    scalar t3 = 1.0 + t1 * alpha;
    scalar t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    out[i] = t4;
    outdx[i] = 0.0;
    outdy[i] = 0.0; // not important
  }
}

//// weak formulation NEWTON  //////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_0_0(int n, double *wt, Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* dodt = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / tau
                      +  vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]
                      - dodt->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_0_0_surf(int n, double *wt, Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (kappa * vj->val[i] * vi->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_0_1(int n, double *wt, Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* dodc = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (- dodc->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_1_0(int n, double *wt, Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* dodt = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( dodt->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_1_1(int n, double *wt, Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* dodc = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / tau
                      +  (vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]) / Le
                      + dodc->val[i] * vj->val[i] * vi->val[i] );
  return result;
}


template<typename Real, typename Scalar>
Scalar newton_linear_form_0(int n, double *wt, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* titer = ext->fn[0];
  Func<Real>* tprev1 = ext->fn[1];
  Func<Real>* tprev2 = ext->fn[2];
  Func<Real>* omega = ext->fn[3];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( (3.0 * titer->val[i] - 4.0 * tprev1->val[i] + tprev2->val[i]) * vi->val[i] / (2.0 * tau) +
                        (titer->dx[i] * vi->dx[i] + titer->dy[i] * vi->dy[i]) -
                        omega->val[i] * vi->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_linear_form_0_surf(int n, double *wt, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (kappa * ext->fn[0]->val[i] * vi->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_linear_form_1(int n, double *wt, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* citer = ext->fn[0];
  Func<Real>* cprev1 = ext->fn[1];
  Func<Real>* cprev2 = ext->fn[2];
  Func<Real>* omega = ext->fn[3];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( (3.0 * citer->val[i] - 4.0 * cprev1->val[i] + cprev2->val[i]) * vi->val[i] / (2.0 * tau) +
                        (citer->dx[i] * vi->dx[i] + citer->dy[i] * vi->dy[i]) / Le +
                        omega->val[i] * vi->val[i]);
  return result;
}

//// main //////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &basemesh);
  basemesh.refine_all_elements();
  basemesh.refine_all_elements();
  mesh.copy(&basemesh);

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // initialize meshes
  H1Space tspace(&mesh, &shapeset);
  H1Space cspace(&mesh, &shapeset);
  tspace.set_bc_types(bc_types);
  tspace.set_bc_values(temp_bc_values);
  cspace.set_bc_types(bc_types);
  tspace.set_uniform_order(init_order);
  cspace.set_uniform_order(init_order);

  // initial conditions
  Solution tprev1, cprev1, tprev2, cprev2, titer, citer, tsln, csln;
  tprev1.set_exact(&mesh, temp_ic);  cprev1.set_exact(&mesh, conc_ic);
  tprev2.set_exact(&mesh, temp_ic);  cprev2.set_exact(&mesh, conc_ic);
  titer.set_exact(&mesh, temp_ic);   citer.set_exact(&mesh, conc_ic);
  DXDYFilter omega(omega_fn, &titer, &citer);
  DXDYFilter omega_dt(omega_dt_fn, &titer, &citer);
  DXDYFilter omega_dc(omega_dc_fn, &titer, &citer);

  ScalarView rview("Reaction rate", 0, 0, 1600, 460);
  OrderView ov("Orders", 0, 500, 1600, 460);

  int ndofs = 0;
  ndofs += tspace.assign_dofs(ndofs);
  ndofs += cspace.assign_dofs(ndofs);

  WeakForm wf(2);
  wf.add_biform(0, 0, callback(newton_bilinear_form_0_0), UNSYM, ANY, 1, &omega_dt);
  wf.add_biform_surf(0, 0, callback(newton_bilinear_form_0_0_surf), 3);
  wf.add_biform(0, 1, callback(newton_bilinear_form_0_1), UNSYM, ANY, 1, &omega_dc);
  wf.add_biform(1, 0, callback(newton_bilinear_form_1_0), UNSYM, ANY, 1, &omega_dt);
  wf.add_biform(1, 1, callback(newton_bilinear_form_1_1), UNSYM, ANY, 1, &omega_dc);
  wf.add_liform(0, callback(newton_linear_form_0), ANY, 4, &titer, &tprev1, &tprev2, &omega);
  wf.add_liform_surf(0, callback(newton_linear_form_0_surf), 3, 1, &titer);
  wf.add_liform(1, callback(newton_linear_form_1), ANY, 4, &citer, &cprev1, &cprev2, &omega);

  UmfpackSolver umfpack;
  NonlinSystem ns(&wf, &umfpack);
  ns.set_spaces(2, &tspace, &cspace);
  ns.set_pss(1, &pss);
  ns.set_ic(&tprev1, &cprev1, &titer, &citer);

  double total_time = 0.0;
  for (int it = 1; total_time <= 60.0; it++)
  {
    info("\n*** Time iteration %d, t = %g s ***", it, total_time + tau);

    int nit = 0;
    bool done = false;
    double nonlin_err;
    do
    {
      info("\n*** Nonlinear iteration %d ***\n", nit++);
      omega.reinit();  omega_dt.reinit();  omega_dc.reinit();

      ns.assemble();
      ns.solve(2, &tsln, &csln);
      nonlin_err = ns.get_residuum_l2_norm();
      info("Nonlinear error: %g\n", nonlin_err);
      if (nonlin_err < 1e-4 || nit > 5) done = true;
      titer = tsln; citer = csln;
    }
    while (!done);

    // visualization
    DXDYFilter omega_view(omega_fn, &titer, &citer);
    rview.set_min_max_range(0.0,2.0);
    rview.show(&omega_view);

    total_time += tau;

    tprev2.copy(&tprev1);
    cprev2.copy(&cprev1);
    tprev1.copy(&titer);
    cprev1.copy(&citer);
  }

  View::wait("Waiting for all views to be closed.");
  return 0;
}
