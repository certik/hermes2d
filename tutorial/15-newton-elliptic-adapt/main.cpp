#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"
#include "function.h"

//  This example shows how the Newton's method can be combined with
//  automatic adaptivity.
//
//  PDE: stationary heat transfer equation with nonlinear thermal
//  conductivity, - div[lambda(u)grad u] = 0
//
//  Domain: unit square (-10,10)^2
//
//  BC: Dirichlet, see function dir_lift() below.
//
//  The following parameters can be changed:

const int P_INIT = 1;                      // Initial polynomial degree
const int PROJ_TYPE = 1;                   // For the projection of the initial condition
                                           // on the initial mesh: 1 = H1 projection, 0 = L2 projection
const int INIT_GLOB_REF_NUM = 1;           // Number of initial uniform mesh refinements
const int INIT_BDY_REF_NUM = 0;            // Number of initial refinements towards boundary

const double THRESHOLD = 0.2;              // This is a quantitative parameter of the adapt(...) function and
                                           // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                    // Adaptive strategy:
                                           // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                           //   error is processed. If more elements have similar errors, refine
                                           //   all to keep the mesh symmetric.
                                           // STRATEGY = 1 ... refine all elements whose error is larger
                                           //   than THRESHOLD times maximum element error.
                                           // STRATEGY = 2 ... refine all elements whose error is larger
                                           //   than THRESHOLD.
                                           // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const RefinementSelectors::AllowedCandidates ADAPT_TYPE = RefinementSelectors::H2DRS_CAND_HP;         // Type of automatic adaptivity.
const bool ISO_ONLY = false;               // Isotropic refinement flag (concerns quadrilateral elements only).
                                           // ISO_ONLY = false ... anisotropic refinement of quad elements
                                           // is allowed (default),
                                           // ISO_ONLY = true ... only isotropic refinements of quad elements
                                           // are allowed.
const int MESH_REGULARITY = -1;            // Maximum allowed level of hanging nodes:
                                           // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                           // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                           // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                           // Note that regular meshes are not supported, this is due to
                                           // their notoriously bad performance.
const double ERR_STOP = 1.0;               // Stopping criterion for adaptivity (rel. error tolerance between the
                                           // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.
const double NEWTON_TOL_COARSE = 1e-6;     // Stopping criterion for the Newton's method on coarse mesh
const double NEWTON_TOL_FINE = 1e-6;       // Stopping criterion for the Newton's method on fine mesh
const int NEWTON_MAX_ITER = 100;           // Maximum allowed number of Newton iterations
const bool NEWTON_ON_COARSE_MESH = false;  // true... Newton is done on coarse mesh in every adaptivity step
                                           // false...Newton is done on coarse mesh only once, then projection
                                           // of the fine mesh solution to coarse mesh is used

// Thermal conductivity (temperature-dependent)
// Note: for any u, this function has to be positive
template<typename Real>
Real lam(Real u)
{
  return 1 + pow(u, 4);
}

// Derivative of the thermal conductivity with respect to 'u'
template<typename Real>
Real dlam_du(Real u) {
  return 4*pow(u, 3);
}

// This function is used to define Dirichlet boundary conditions
double dir_lift(double x, double y, double& dx, double& dy) {
  dx = (y+10)/10.;
  dy = (x+10)/10.;
  return (x+10)*(y+10)/100.;
}

// This function will be projected on the initial mesh and
// used as initial guess for the Newton's method
scalar init_guess(double x, double y, double& dx, double& dy)
{
  // using the Dirichlet lift elevated by two
  double val = dir_lift(x, y, dx, dy) + 2;
  return val;
}

// Boundary condition type (essential = Dirichlet)
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Dirichlet boundary condition values
scalar bc_values(int marker, double x, double y)
{
  double dx, dy;
  return dir_lift(x, y, dx, dy);
}

// Heat sources (can be a general function of 'x' and 'y')
template<typename Real>
Real heat_src(Real x, Real y)
{
  return 1.0;
}

// Jacobian matrix
template<typename Real, typename Scalar>
Scalar jac(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (dlam_du(u_prev->val[i]) * u->val[i] * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i])
                       + lam(u_prev->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));

  return result;
}

// Fesidual vector
template<typename Real, typename Scalar>
Scalar res(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (lam(u_prev->val[i]) * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i])
		       - heat_src(e->x[i], e->y[i]) * v->val[i]);
  return result;
}

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // initial mesh refinements
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // enumerate degrees of freedom
  int ndof = assign_dofs(&space);

  // solutions for the Newton's iteration and adaptivity
  Solution u_prev, sln_coarse, sln_fine;

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(jac), H2D_UNSYM, H2D_ANY, 1, &u_prev);
  wf.add_liform(0, callback(res), H2D_ANY, 1, &u_prev);

  // initialize the nonlinear system and solver
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof, graph_cpu;

  // time measurement
  TimePeriod cpu_time;

  // project the function init_guess() on the coarse mesh
  // to obtain initial guess u_prev for the Newton's method
  nls.set_ic(init_guess, &mesh, &u_prev, PROJ_TYPE);

  // visualisation
  ScalarView sview_coarse("Coarse mesh solution", 0, 0, 500, 400); // coarse mesh solution
  OrderView oview_coarse("Coarse mesh", 520, 0, 450, 400);         // coarse mesh
  ScalarView sview_fine("Fine mesh solution", 990, 0, 500, 400);   // fine mesh solution
  OrderView oview_fine("Fine mesh", 1510, 0, 450, 400);            // fine mesh

  // Newton's loop on the coarse mesh
  info("---- Solving on coarse mesh:\n");
  if (!nls.solve_newton_1(&u_prev, NEWTON_TOL_COARSE, NEWTON_MAX_ITER)) error("Newton's method did not converge.");

  // store the result in sln_coarse
  sln_coarse.copy(&u_prev);

  // prepare selector
  RefinementSelectors::H1NonUniformHP selector(ISO_ONLY, ADAPT_TYPE, 1.0, H2DRS_DEFAULT_ORDER, &shapeset);

  // adaptivity loop
  bool done = false;
  int a_step = 0;
  do {
    double err_est;

    a_step++;
    info("!---- Adaptivity step %d ---------------------------------------------", a_step);

    // show coarse mesh and solution
    sview_coarse.show(&sln_coarse);
    oview_coarse.show(&space);

    info("---- Solving on fine mesh ---------------------------------------------");

    // do not measure showing of results
    cpu_time.tick(H2D_SKIP);

    // Setting initial guess for the Newton's method on the fine mesh
    RefNonlinSystem rnls(&nls);
    rnls.prepare();
    if (a_step == 1) rnls.set_ic(&sln_coarse, &u_prev, PROJ_TYPE);
    else rnls.set_ic(&sln_fine, &u_prev, PROJ_TYPE);

    // Newton's loop on the fine mesh
    if (!rnls.solve_newton_1(&u_prev, NEWTON_TOL_FINE, NEWTON_MAX_ITER)) error("Newton's method did not converge.");

    // stote the fine mesh solution in sln_fine
    sln_fine.copy(&u_prev);

    // time measurement
    cpu_time.tick();

    // show fine mesh and solution
    sview_fine.show(&sln_fine);
    oview_fine.show(rnls.get_space(0));

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // calculate element errors and total error estimate
    H1AdaptHP hp(1, &space);
    err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;

    // time measurement
    cpu_time.tick();

    // report results
    info("Error estimate: %g%%", err_est);

    // add entry to DOF convergence graph
    graph_dof.add_values(space.get_num_dofs(), err_est);
    graph_dof.save("conv_dof.dat");

    // add entry to CPU convergence graph
    graph_cpu.add_values(cpu_time.accumulated(), err_est);
    graph_cpu.save("conv_cpu.dat");

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(THRESHOLD, STRATEGY, &selector, MESH_REGULARITY);
      ndof = assign_dofs(&space);
      if (ndof >= NDOF_STOP) done = true;

      // project the fine mesh solution on the new coarse mesh
      info("---- Projecting fine mesh solution on new coarse mesh -----------------");
      nls.set_ic(&sln_fine, &u_prev, PROJ_TYPE);

      if (NEWTON_ON_COARSE_MESH) {
        // Newton's loop on the coarse mesh
        info("---- Solving on coarse mesh -------------------------------------------");
        if (!nls.solve_newton_1(&u_prev, NEWTON_TOL_COARSE, NEWTON_MAX_ITER)) error("Newton's method did not converge.");
      }

      // store the result in sln_coarse
      sln_coarse.copy(&u_prev);
    }

    // time measurement
    cpu_time.tick();
  }
  while (!done);
  cpu_time.tick();
  verbose("Total running time: %g s", cpu_time.accumulated());

  // wait for keyboard or mouse input
  View::wait();
  return 0;
}

