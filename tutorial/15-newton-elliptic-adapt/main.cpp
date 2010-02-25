#include "hermes2d.h"
#include "solver_umfpack.h"
#include "function.h"

//  This example shows how to combine the Newton's method with 
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

const int P_INIT = 1;                // Initial polynomial degree
const int PROJ_TYPE = 1;             // For the projection of the initial condition 
                                     // on the initial mesh: 1 = H1 projection, 0 = L2 projection
const int INIT_GLOB_REF_NUM = 1;     // Number of initial uniform mesh refinements
const int INIT_BDY_REF_NUM = 0;      // Number of initial refinements towards boundary

const double THRESHOLD = 0.2;        // This is a quantitative parameter of the adapt(...) function and
                                     // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;              // Adaptive strategy:
                                     // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                     //   error is processed. If more elements have similar errors, refine
                                     //   all to keep the mesh symmetric.
                                     // STRATEGY = 1 ... refine all elements whose error is larger
                                     //   than THRESHOLD times maximum element error.
                                     // STRATEGY = 2 ... refine all elements whose error is larger
                                     //   than THRESHOLD.
                                     // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int ADAPT_TYPE = 0;            // Type of automatic adaptivity:
                                     // ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                     // ADAPT_TYPE = 1 ... adaptive h-FEM,
                                     // ADAPT_TYPE = 2 ... adaptive p-FEM.
const bool ISO_ONLY = false;         // Isotropic refinement flag (concerns quadrilateral elements only).
                                     // ISO_ONLY = false ... anisotropic refinement of quad elements
                                     // is allowed (default),
                                     // ISO_ONLY = true ... only isotropic refinements of quad elements
                                     // are allowed.
const int MESH_REGULARITY = -1;      // Maximum allowed level of hanging nodes:
                                     // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                     // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                     // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                     // Note that regular meshes are not supported, this is due to
                                     // their notoriously bad performance.
const double ERR_STOP = 0.001;       // Stopping criterion for adaptivity (rel. error tolerance between the
                                     // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;         // Adaptivity process stops when the number of degrees of freedom grows
                                     // over this limit. This is to prevent h-adaptivity to go on forever.
const double NEWTON_TOL = 1e-6;      // Stopping criterion for the Newton's method on coarse mesh
const double NEWTON_TOL_REF = 1e-6;  // Stopping criterion for the Newton's method on fine mesh

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
scalar init_cond(double x, double y, double& dx, double& dy)
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
  mesh.refine_towards_boundary(1,INIT_BDY_REF_NUM);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
  space.assign_dofs();

  // previous solution for the Newton's iteration
  Solution u_prev;

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(jac), UNSYM, ANY, 1, &u_prev);
  wf.add_liform(0, callback(res), ANY, 1, &u_prev);

  // initialize the nonlinear system and solver
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof, graph_cpu;

  // project the function init_cond() on the mesh 
  // to obtain initial guess u_prev for the Newton's method
  nls.set_ic(init_cond, &mesh, &u_prev, PROJ_TYPE);

  // visualise the initial ocndition
  ScalarView view("Initial condition", 0, 0, 700, 600);
  view.show(&u_prev);
  OrderView oview("Initial mesh", 720, 0, 700, 600);
  oview.show(&space);
  //printf("Click into the image window and press any key to proceed.\n");
  //view.wait_for_keypress();

  // adaptivity loop
  double cpu = 0.0, err_est;
  int a_step = 1;
  bool done = false;
  do {

    a_step++;

    // Newton's loop on the coarse mesh
    int it = 1;
    double res_l2_norm;
    Solution sln_coarse;
    do
    {
      info("\n---- Adapt step %d, Newton iter %d (coarse mesh) ---------------------------------\n", a_step, it++);
      printf("ndof = %d\n", space.get_num_dofs());

      // time measurement
      begin_time();

      // assemble the Jacobian matrix and residual vector, 
      // solve the system
      nls.assemble();
      nls.solve(1, &sln_coarse);

      // calculate the l2-norm of residual vector
      res_l2_norm = nls.get_residuum_l2_norm();
      info("Residuum L2 norm: %g\n", res_l2_norm);

      // time measurement
      cpu += end_time();

      // visualise the solution
      char title[100];
      sprintf(title, "Temperature (coarse mesh), Newton iteration %d", it-1);
      view.set_title(title);
      view.show(&sln_coarse);
      sprintf(title, "Coarse mesh, Newton iteration %d", it-1);
      oview.set_title(title);
      oview.show(&space);
      //printf("Click into the image window and press any key to proceed.\n");
      //view.wait_for_keypress();

      // save the new solution as "previous" for the 
      // next Newton's iteration
      u_prev.copy(&sln_coarse);
    }
    while (res_l2_norm > NEWTON_TOL);

    // Setting initial guess for the Newton's method on the fine mesh
    Solution sln_fine, u_prev_fine;
    RefNonlinSystem rs(&nls);
    rs.prepare();
    rs.set_ic(&u_prev, &u_prev);

    // Newton's loop on the fine mesh
    it = 1;
    do {
      info("\n---- Adapt step %d, Newton iter %d (fine mesh) ---------------------------------\n", a_step, it++);

      // time measurement
      begin_time();

      // assemble the Jacobian matrix and residual vector, 
      // solve the system
      rs.assemble();
      rs.solve(1, &sln_fine);

      // calculate the l2-norm of residual vector
      res_l2_norm = rs.get_residuum_l2_norm();
      info("Residuum L2 norm: %g\n", res_l2_norm);

      // time measurement
      cpu += end_time();

      // visualise the solution
      char title[100];
      sprintf(title, "Temperature (fine mesh), Newton iteration %d", it-1);
      view.set_title(title);
      view.show(&sln_fine);
      sprintf(title, "Fine mesh, Newton iteration %d", it-1);
      oview.set_title(title);
      oview.show(rs.get_ref_space(0));
      //printf("Click into the image window and press any key to proceed.\n");
      //view.wait_for_keypress();

      u_prev.copy(&sln_fine);
    } while (res_l2_norm > NEWTON_TOL_REF);

    // time measurement
    begin_time();

    // calculate element errors and total error estimate
    H1OrthoHP hp(1, &space);
    err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;
    info("Error estimate: %g%%", err_est);

    // add entry to DOF convergence graph
    graph_dof.add_values(space.get_num_dofs(), err_est);
    graph_dof.save("conv_dof.dat");

    // add entry to CPU convergence graph
    graph_cpu.add_values(cpu, err_est);
    graph_cpu.save("conv_cpu.dat");

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE, ISO_ONLY, MESH_REGULARITY);
      int ndof = space.assign_dofs();
      if (ndof >= NDOF_STOP) done = true;
    }

    // time measurement
    cpu += end_time();
  }
  while (!done);
  verbose("Total running time: %g sec", cpu);

  // wait for keyboard or mouse input
  View::wait();
  return 0;
}

