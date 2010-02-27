#include "hermes2d.h"
#include "solver_umfpack.h"

//  This example is a very simple flame propagation model (laminar flame,
//  zero flow velocity), and its purpose is to show how the Newton's method
//  is applied to a time-dependent two-equation system.
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

const int P_INIT = 2;                  // Initial polynomial degree
const double TAU = 0.5;                // Time step
const double NEWTON_TOL = 1e-4;        // Stopping criterion for the Newton's method on coarse mesh
const int NEWTON_MAX_ITER = 10;        // Maximum allowed number of Newton iterations

// Problem constants
const double Le    = 1.0;
const double alpha = 0.8;
const double beta  = 10.0;
const double kappa = 0.1;
const double x1    = 9.0;

// Boundary and initial conditions
int bc_types(int marker)
  { return (marker == 1) ? BC_ESSENTIAL : BC_NATURAL; }

scalar temp_bc_values(int marker, double x, double y)
  { return (marker == 1) ? 1.0 : 0; }

scalar temp_ic(double x, double y, scalar& dx, scalar& dy)
  { return (x <= x1) ? 1.0 : exp(x1 - x); }

scalar conc_ic(double x, double y, scalar& dx, scalar& dy)
  { return (x <= x1) ? 0.0 : 1.0 - exp(Le*(x1 - x)); }

// Weak forms, definition of reaction rate omega
# include "forms.cpp"

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // initial mesh refinements
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create H1 spaces
  H1Space tspace(&mesh, &shapeset);
  H1Space cspace(&mesh, &shapeset);
  tspace.set_bc_types(bc_types);
  tspace.set_bc_values(temp_bc_values);
  cspace.set_bc_types(bc_types);
  tspace.set_uniform_order(P_INIT);
  cspace.set_uniform_order(P_INIT);
  int ndofs = 0;
  ndofs += tspace.assign_dofs(ndofs);
  ndofs += cspace.assign_dofs(ndofs);

  // solutions for the Newton's iteration and time stepping
  Solution tprev1, cprev1, tprev2, cprev2, titer, citer, tsln, csln;

  // setting initial conditions
  tprev1.set_exact(&mesh, temp_ic);  cprev1.set_exact(&mesh, conc_ic);
  tprev2.set_exact(&mesh, temp_ic);  cprev2.set_exact(&mesh, conc_ic);
  titer.set_exact(&mesh, temp_ic);   citer.set_exact(&mesh, conc_ic);

  // defining filters for the reaction rate omega
  DXDYFilter omega(omega_fn, &titer, &citer);
  DXDYFilter omega_dt(omega_dt_fn, &titer, &citer);
  DXDYFilter omega_dc(omega_dc_fn, &titer, &citer);

  // visualization
  ScalarView rview("Reaction rate", 0, 0, 1600, 460);

  // initialize the weak formulation
  WeakForm wf(2);
  wf.add_biform(0, 0, callback(newton_bilinear_form_0_0), UNSYM, ANY, 1, &omega_dt);
  wf.add_biform_surf(0, 0, callback(newton_bilinear_form_0_0_surf), 3);
  wf.add_biform(0, 1, callback(newton_bilinear_form_0_1), UNSYM, ANY, 1, &omega_dc);
  wf.add_biform(1, 0, callback(newton_bilinear_form_1_0), UNSYM, ANY, 1, &omega_dt);
  wf.add_biform(1, 1, callback(newton_bilinear_form_1_1), UNSYM, ANY, 1, &omega_dc);
  wf.add_liform(0, callback(newton_linear_form_0), ANY, 4, &titer, &tprev1, &tprev2, &omega);
  wf.add_liform_surf(0, callback(newton_linear_form_0_surf), 3, 1, &titer);
  wf.add_liform(1, callback(newton_linear_form_1), ANY, 4, &citer, &cprev1, &cprev2, &omega);

  // initialize the nonlinear system and solver
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(2, &tspace, &cspace);
  nls.set_pss(1, &pss);
  nls.set_ic(&tprev1, &cprev1, &titer, &citer);

  // time stepping loop
  double total_time = 0.0;
  for (int it = 1; total_time <= 60.0; it++)
  {
    info("\n**** Time step %d, t = %g s:\n", it, total_time);

    nls.solve_newton_2(&titer, &citer, NEWTON_TOL, NEWTON_MAX_ITER, 
                       &omega, &omega_dt, &omega_dc);

    // visualization
    DXDYFilter omega_view(omega_fn, &titer, &citer);
    rview.set_min_max_range(0.0,2.0);
    char title[100];
    sprintf(title, "Reaction rate, t = %g", total_time);
    rview.set_title(title);
    rview.show(&omega_view);

    total_time += TAU;

    tprev2.copy(&tprev1);
    cprev2.copy(&cprev1);
    tprev1.copy(&titer);
    cprev1.copy(&citer);
  }

  View::wait("Waiting for all views to be closed.");
  return 0;
}
