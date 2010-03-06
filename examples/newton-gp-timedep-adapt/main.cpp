#include "hermes2d.h"
#include "solver_umfpack.h"

//  This example shows how to combine automatic adaptivity with the Newton's
//  method for a nonlinear complex-valued time-dependent PDE (the Gross-Pitaevski
//  equation describing the behavior of Einstein-Bose quantum gases)
//  discretized implicitly in time (via implicit Euler or Crank-Nicolson).
//
//  PDE: non-stationary complex Gross-Pitaevski equation
//  describing resonances in Bose-Einstein condensates
//
//  ih \partial \psi/\partial t = -h^2/(2m) \Delta \psi +
//  g \psi |\psi|^2 + 1/2 m \omega^2 (x^2 + y^2) \psi
//
//  square (-1, 1)^2
//
//  BC:  homogeneous Dirichlet everywhere on the boundary
//
//  Time-stepping: either implicit Euler or Crank-Nicolson

const int P_INIT = 1;                    // Initial polynomial degree
const int PROJ_TYPE = 1;                 // For the projection of the initial condition 
                                         // on the initial mesh: 1 = H1 projection, 0 = L2 projection
const int INIT_REF_NUM = 2;              // Number of initial uniform refinements
const int TIME_DISCR = 2;                // 1 for implicit Euler, 2 for Crank-Nicolson
const double T_FINAL = 200.0;            // Time interval length
const double TAU = 0.01;                 // Time step

// Newton's method
const double NEWTON_TOL_COARSE = 0.01;   // Stopping criterion for Newton on coarse mesh
const double NEWTON_TOL_FINE = 0.05;     // Stopping criterion for Newton on fine mesh
                                         // (the ref. solution does not have to be super-accurate)
const int NEWTON_MAX_ITER = 100;         // Maximum allowed number of Newton iterations

// Adaptivity
const int UNREF_FREQ = 1;                // Every UNREF_FREQ time step the mesh is unrefined
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int ADAPT_TYPE = 0;                // Type of automatic adaptivity:
                                         // ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                         // ADAPT_TYPE = 1 ... adaptive h-FEM,
                                         // ADAPT_TYPE = 2 ... adaptive p-FEM.
const bool ISO_ONLY = false;             // Isotropic refinement flag (concerns quadrilateral elements only).
                                         // ISO_ONLY = false ... anisotropic refinement of quad elements
                                         // is allowed (default),
                                         // ISO_ONLY = true ... only isotropic refinements of quad elements
                                         // are allowed.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of 
                                         // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int MAX_P = 5;                     // Maximum polynomial order allowed in hp-adaptivity
                                         // had to be limited due to complicated integrals
const double ERR_STOP = 1.0;             // Stopping criterion for hp-adaptivity
                                         // (relative error between reference and coarse solution in percent)
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows
                                         // over this limit. This is to prevent h-adaptivity to go on forever.
const int SHOW_MESHES_IN_TIME_STEP = 1;  // If nonzero, all meshes during every time step are
                                         // shown, else only meshes at the end of every time step.

// Problem constants
const double H = 1;                      // Planck constant 6.626068e-34;
const double M = 1;                      // mass of boson
const double G = 1;                      // coupling constant
const double OMEGA = 1;                  // frequency


// Initial conditions
scalar fn_init(double x, double y, scalar& dx, scalar& dy)
{
  scalar val = exp(-20*(x*x + y*y));
  dx = val * (-40.0*x);
  dy = val * (-40.0*y);
  return val;
}

// Boundary condition types
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Boundary condition values
scalar bc_values(int marker, double x, double y)
{
 return 0;
}

// Weak forms
# include "integrals.cpp"

// Implicit Euler method (1st-order in time)
template<typename Real, typename Scalar>
Scalar residuum_euler(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{  return F_euler(n, wt, v, e, ext);  }
template<typename Real, typename Scalar>
Scalar jacobian_euler(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{  return J_euler(n, wt, u, v, e, ext);  }

// Crank-Nicolson (2nd-order in time)
template<typename Real, typename Scalar>
Scalar residuum_cranic(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{  return F_cranic(n, wt, v, e, ext);  }
template<typename Real, typename Scalar>
Scalar jacobian_cranic(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{  return J_cranic(n, wt, u, v, e, ext);  }

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("square.mesh", &basemesh);

  // initial mesh refinements
  for(int i = 0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
  space.assign_dofs();

  // solutions for the Newton's iteration and adaptivity
  Solution Psi_prev_time, Psi_prev_newton, sln_coarse, sln_fine; 

  // initialize the weak formulation
  WeakForm wf(1);
  if(TIME_DISCR == 1) {
    wf.add_biform(0, 0, callback(jacobian_euler), UNSYM, ANY, 1, &Psi_prev_newton);
    wf.add_liform(0, callback(residuum_euler), ANY, 2, &Psi_prev_newton, &Psi_prev_time);
  }
  else {
    wf.add_biform(0, 0, callback(jacobian_cranic), UNSYM, ANY, 1, &Psi_prev_newton);
    wf.add_liform(0, callback(residuum_cranic), ANY, 2, &Psi_prev_newton, &Psi_prev_time);
  }

  // initialize the nonlinear system and solver
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  // visualization
  char title[100];
  ScalarView view("", 0, 0, 600, 500);
  view.fix_scale_width(80);
  ScalarView magview("", 0, 0, 600, 500);
  magview.fix_scale_width(80);
  OrderView ordview("", 610, 0, 600, 500);
  ordview.fix_scale_width(80);

  // DOF and CPU convergence graphs
  SimpleGraph graph_time_dof, graph_time_err;

  // set initial condition at zero time level
  Psi_prev_time.set_exact(&mesh, fn_init);
  Psi_prev_newton.set_exact(&mesh, fn_init);
  nls.set_ic(&Psi_prev_newton, &Psi_prev_newton, PROJ_TYPE);
  sln_coarse.copy(&Psi_prev_time);

  // time stepping loop
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int n = 1; n <= nstep; n++)
  {
    info("\n---- Time step %d:\n", n);

    // periodical global derefinements
    if (n % UNREF_FREQ == 0) {
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);
      space.assign_dofs();
    }

    // adaptivity loop
    int a_step = 0;
    bool done = false;
    double err_est;
    do
    {
      info("---- Time step %d, adaptivity step %d, coarse mesh solution:\n", n, ++a_step);

      // Newton's method on coarse mesh
      if (!nls.solve_newton_1(&Psi_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER)) error("Newton's method did not converge.");
      sln_coarse.copy(&Psi_prev_newton);

      info("---- Time step %d, adaptivity step %d, fine mesh solution:\n", n, ++a_step);

      // reference nonlinear system
      RefNonlinSystem rnls(&nls);
      rnls.prepare();

      // set initial condition for the Newton's method on the fine mesh
      if (n == 1 || a_step == 1) rnls.set_ic(&sln_coarse, &Psi_prev_newton);
      else rnls.set_ic(&sln_fine, &Psi_prev_newton);

      // Newton's method on fine mesh
      if (!rnls.solve_newton_1(&Psi_prev_newton, NEWTON_TOL_FINE, NEWTON_MAX_ITER)) error("Newton's method did not converge.");
      sln_fine.copy(&Psi_prev_newton);

      // visualization of intermediate solution
      // and mesh during adaptivity
      if(SHOW_MESHES_IN_TIME_STEP) {
        sprintf(title, "Solution magnitude, time level %d", n);
        magview.set_title(title);
        AbsFilter mag(&Psi_prev_newton);
        magview.show(&mag);            
        sprintf(title, "hp-mesh, time level %d", n);
        ordview.set_title(title);
        ordview.show(&space);          
      }

      // calculate element errors and total error estimate
      H1OrthoHP hp(1, &space);
      err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;   // relative h1-error in percent
      info("Error estimate: %g%%", err_est);

      // if err_est too large, adapt the mesh
      if (err_est < ERR_STOP) done = true;
      else {
        hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE, ISO_ONLY, MESH_REGULARITY, CONV_EXP, MAX_P);
        int ndof = space.assign_dofs();
        if (ndof >= NDOF_STOP) done = true;

        // update initial guess Psi_prev_newton for the Newton's method
        // on the new coarse mesh
        nls.set_ic(&Psi_prev_newton, &Psi_prev_newton, PROJ_TYPE);
      }
    }
    while (!done);

    // show magnitude of the solution
    sprintf(title, "Magnitude, time level %d", n);
    magview.set_title(title);
    AbsFilter mag(&sln_fine);
    magview.show(&mag);      

    // show hp-mesh
    sprintf(title, "hp-mesh, time level %d", n);
    ordview.set_title(title);
    ordview.show(&space);

    // add entries to convergence graphs
    graph_time_err.add_values(n*TAU, err_est);
    graph_time_err.save("time_error.dat");
    graph_time_dof.add_values(n*TAU, space.get_num_dofs());
    graph_time_dof.save("time_dof.dat");

    // copy result of the Newton's iteration into Psi_prev_time
    Psi_prev_time.copy(&Psi_prev_newton);
  }

  // wait for keyboard or mouse input
  View::wait();
  return 0;
}
