#include "hermes2d.h"
#include "solver_umfpack.h"

//  This example shows how to combine automatic adaptivity with the Newton's
//  method for a nonlinear time-dependent PDE discretized implicitly in time
//  (using implicit Euler or Crank-Nicolson). Unrefinements are allowed.
//  Some problem parameters can be changed below.
//
//  PDE: non-stationary heat transfer with nonlinear thermal conductivity
//  HEATCAP*dT/dt - div[lambda(T)grad T] = 0
//
//  hp-adaptivity with dynamical meshes
//
//  Domain: square
//
//  BC:  T = 100 on the left, top and bottom edges
//       dT/dn = 0 on the right edge
//
//  The following parameters can be changed:

// adaptivity parameters
const int P_INIT = 2;             // Initial polynomial degree of all mesh elements.
const int REF_INIT = 2;           // Number of initial refinements.
const int UNREF_FREQ = 1;         // Every UNREF_FREQth time step the mesh is unrefined.
const double SPACE_H1_TOL = 0.5;  // Stopping criterion for hp-adaptivity
                                  // (relative error between reference and coarse solution in percent).
const double THRESHOLD = 0.3;     // This is a quantitative parameter of the adapt(...) function and
                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;           // Adaptive strategy:
                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                  //   error is processed. If more elements have similar errors, refine
                                  //   all to keep the mesh symmetric.
                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                  //   than THRESHOLD times maximum element error.
                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                  //   than THRESHOLD.
                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int ADAPT_TYPE = 0;         // Type of automatic adaptivity:
                                  // ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                  // ADAPT_TYPE = 1 ... adaptive h-FEM,
                                  // ADAPT_TYPE = 2 ... adaptive p-FEM.
const bool ISO_ONLY = false;      // Isotropic refinement flag (concerns quadrilateral elements only).
                                  // ISO_ONLY = false ... anisotropic refinement of quad elements
                                  // is allowed (default),
                                  // ISO_ONLY = true ... only isotropic refinements of quad elements
                                  // are allowed.
const int MESH_REGULARITY = -1;   // Maximum allowed level of hanging nodes:
                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                  // Note that regular meshes are not supported, this is due to
                                  // their notoriously bad performance.
const double ERR_STOP = 0.15;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 50000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// time discretization parameters
const int TIME_DISCR = 1;         // 1 for implicit Euler, 2 for Crank-Nicolson
const int PROJ_TYPE = 1;          // 1 for H1 projections, 0 for L2 projections
const double HEATCAP = 1e6;       // heat capacity
const double TAU = 5;             // time step
const double T_FINAL = 600;       // time interval length

// Newton parameters
const double NEWTON_TOL_COARSE = 0.05;   // stopping criterion for Newton on coarse mesh
const double NEWTON_TOL_REF = 0.5;       // stopping criterion for Newton on fine mesh
                                   // (the ref. solution does not to be super accurate)

// thermal conductivity (nonlinear, temperature-dependent
// for any u, this function has to be positive in the entire
// to ensure solvability!
template<typename Real>
Real lam(Real T) { return 10 + 0.1*pow(T, 2); }
template<typename Real>
Real dlam_dT(Real T) { return 0.1*2*pow(T, 1); }

// boundary conditions
int bc_types(int marker)
{
  if (marker == 1) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

scalar bc_values(int marker, double x, double y)
{
  return 100;
}

// Jacobian matrices and residual vectors

// Residuum for the implicit Euler time discretization
template<typename Real, typename Scalar>
Scalar F_euler(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* titer = ext->fn[0];
  Func<Scalar>* tprev = ext->fn[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * (HEATCAP*(titer->val[i] - tprev->val[i]) * v->val[i] / TAU +
                       lam(titer->val[i]) * (titer->dx[i] * v->dx[i] + titer->dy[i] * v->dy[i]));
  return result;
}

template<typename Real, typename Scalar>
Scalar J_euler(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* titer = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (HEATCAP * u->val[i] * v->val[i] / TAU +
                       dlam_dT(titer->val[i]) * u->val[i] * (titer->dx[i] * v->dx[i] + titer->dy[i] * v->dy[i]) +
                       lam(titer->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
  return result;
}

// Residuum for the implicit Euler time discretization
template<typename Real, typename Scalar>
Scalar F_cranic(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* titer = ext->fn[0];
  Func<Scalar>* tprev = ext->fn[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * (HEATCAP * (titer->val[i] - tprev->val[i]) * v->val[i] / TAU +
                       0.5 * lam(titer->val[i]) * (titer->dx[i] * v->dx[i] + titer->dy[i] * v->dy[i]) +
                       0.5 * lam(tprev->val[i]) * (tprev->dx[i] * v->dx[i] + tprev->dy[i] * v->dy[i]));
  return result;
}

template<typename Real, typename Scalar>
Scalar J_cranic(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* titer = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (HEATCAP * u->val[i] * v->val[i] / TAU +
                       0.5 * dlam_dT(titer->val[i]) * u->val[i] * (titer->dx[i] * v->dx[i] + titer->dy[i] * v->dy[i]) +
                       0.5 * lam(titer->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
  return result;
}

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("square.mesh", &basemesh);
  for(int i = 0; i < REF_INIT; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);
  mesh.refine_towards_boundary(1,3);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create finite element space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
  space.assign_dofs();

  // enumerate basis functions
  space.assign_dofs();

  Solution Tprev, // previous time step solution, for the time integration method
           Titer; // solution converging during the Newton's iteration

  // initialize the weak formulation
  WeakForm wf(1);
  if(TIME_DISCR == 1) {
    wf.add_biform(0, 0, callback(J_euler), UNSYM, ANY, 1, &Titer);
    wf.add_liform(0, callback(F_euler), ANY, 2, &Titer, &Tprev);
  }
  else {
    wf.add_biform(0, 0, callback(J_cranic), UNSYM, ANY, 1, &Titer);
    wf.add_liform(0, callback(F_cranic), ANY, 2, &Titer, &Tprev);
  }

  // matrix solver
  UmfpackSolver solver;

  // nonlinear system class
  NonlinSystem nls(&wf, &solver);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  // visualize solution and mesh
  ScalarView view("", 0, 0, 700, 600);
  view.fix_scale_width(80);
  OrderView ordview("", 700, 0, 700, 600);

  // error estimate as a function of physical time
  GnuplotGraph graph_err;
  graph_err.set_captions("","Time step","Error");
  graph_err.add_row();

  // error estimate as a function of DOF
  GnuplotGraph graph_dofs;
  graph_dofs.set_captions("","Time step","DOFs");
  graph_dofs.add_row();

  // initial condition at zero time level
  //Tprev.set_const(&mesh, 0.0);
  Tprev.set_dirichlet_lift(&space, &pss);
  Titer.set_dirichlet_lift(&space, &pss);
  nls.set_ic(&Titer, &Titer, PROJ_TYPE);

  // view initial guess for Newton's method
  // satisfies BC conditions
  char title[100];
  sprintf(title, "Initial iteration");
  view.set_title(title);
  view.show(&Titer);
  ordview.show(&space);
  //view.wait_for_keypress(); // this may cause graphics problems

  // time stepping loop
  int nstep = (int)(T_FINAL/TAU + 0.5);
  double cpu = 0.0;
  Solution sln_coarse, sln_fine;
  for(int n = 1; n <= nstep; n++)
  {

    info("\n---- Time step %d -----------------------------------------------------------------", n);

    // time measurement
    begin_time();

    // perform periodic unrefinements
    if (n % UNREF_FREQ == 0) {
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);
      space.assign_dofs();
    }

    // adaptivity loop
    int at = 0, ndofs;
    bool done = false;
    double err_est, cpu;
    do
    {
     info("\n---- Time step %d, adaptivity step %d ---------------------------------------------\n", n, ++at);

      // Newton's loop for coarse mesh solution
      int it = 1;
      double res_l2_norm;
      if (n > 1 || at > 1) nls.set_ic(&sln_fine, &Titer);
      else nls.set_ic(&Titer, &Titer);
      do
      {
        info("\n---- Time step %d, adaptivity step %d, Newton step %d (Coarse mesh solution)-------\n", n, at, it++);

        nls.assemble();
        nls.solve(1, &sln_coarse);

        res_l2_norm = nls.get_residuum_l2_norm();
        info("Residuum L2 norm: %g", res_l2_norm);

        Titer.copy(&sln_coarse);
      }
      while (res_l2_norm > NEWTON_TOL_COARSE);

      // Newton's loop for fine mesh solution
      it = 1;
      RefNonlinSystem rs(&nls);
      rs.prepare();
      if (n > 1 || at > 1) rs.set_ic(&sln_fine, &Titer);
      else rs.set_ic(&Titer, &Titer);
      do
      {
        info("\n---- Time step %d, adaptivity step %d, Newton step %d (Fine mesh solution) --------\n", n, at, it++);

        rs.assemble();
        rs.solve(1, &sln_fine);

        res_l2_norm = rs.get_residuum_l2_norm();
        info("Residuum L2 norm: %g", res_l2_norm);

        Titer.copy(&sln_fine);
      }
      while (res_l2_norm > NEWTON_TOL_REF);

      // calculate error estimate wrt. fine mesh solution
      H1OrthoHP hp(1, &space);
      err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;
      info("Error estimate: %g%", err_est);

      // visualization of solution on the n-th time level
      sprintf(title, "Temperature, time level %d", n);
      //view.set_min_max_range(0,100);
      view.set_title(title);
      //view.show(&Titer);    // to see reference solution
      view.show(&sln_fine);        // to see the solution

      // visualization of mesh on the n-th time level
      sprintf(title, "hp-mesh, time level %d", n);
      ordview.set_title(title);
      ordview.show(&space);   // to see hp-mesh
      //view.wait_for_keypress();

      // if err_est too large, adapt the mesh
      if (err_est < SPACE_H1_TOL) done = true;
      else {
        hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE, ISO_ONLY, MESH_REGULARITY);
        ndofs = space.assign_dofs();
        if (ndofs >= NDOF_STOP) done = true;
      }

      // time measurement
      cpu += end_time();
    }
    while (!done);

    // add entry to both time and DOF error graphs
    graph_err.add_values(0, n, err_est);
    graph_err.save("error.txt");
    graph_dofs.add_values(0, n, space.get_num_dofs());
    graph_dofs.save("dofs.txt");

    // copying result of the Newton's iteration into Tprev
    Tprev.copy(&Titer);
  }

  // time measurement
  cpu += end_time();
  verbose("Total running time: %g sec", cpu);

  // wait for keyboard or mouse input
  View::wait("Waiting for all views to be closed.");
  return 0;
}
