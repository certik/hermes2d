#include "hermes2d.h"
#include "solver_umfpack.h"

//
//  This example demonstrates the employment of the Newton's method to
//  a nonlinear time-dependent PDE discretized implicitly in time
//  (via implicit Euler or Crank-Nicolson). Some problem parameters can
//  be changed below.
//
//  PDE: non-stationary heat transfer with nonlinear thermal conductivity
//  HEATCAP*dT/dt - div[lambda(T)grad T] = 0
//
//  Domain: square
//
//  BC:  T = 100 on the left, top and bottom edges
//       dT/dn = 0 on the right edge
//
//  Time-stepping: either implicit Euler or Crank-Nicolson
//

/********** Problem parameters ***********/

int TIME_DISCR = 2;        // 1 for implicit Euler, 2 for Crank-Nicolson
int PROJ_TYPE = 1;         // 1 for H1 projections, 0 for L2 projections
double HEATCAP = 1e6;      // heat capacity
double TAU = 0.5;          // time step
int NSTEP = 1000;          // number of time steps to do
double NEWTON_TOL = 1e-3;  // convergence criterion for the Newton's method

// thermal conductivity (temperature-dependent
// for any u, this function has to be  positive in the entire domain!
template<typename Real>
Real lam(Real T) { return 10 + 0.1*pow(T, 2); }
template<typename Real>
Real dlam_dT(Real T) { return 0.1*2*pow(T, 1); }

/********** Definition of boundary conditions ***********/

int bc_types(int marker)
{
 if (marker == 4 || marker == 1 || marker == 3) return BC_ESSENTIAL;
//   if (marker == 4) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

scalar bc_values(int marker, double x, double y)
{
 return 100;
// return -4.0 * sqr(y) + 4.0 * y;
}


/********** Definition of Jacobian matrices and residual vectors ***********/

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

// *************************************************************

int main(int argc, char* argv[])
{
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);
  for(int i = 0; i < 5; i++) mesh.refine_all_elements();

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(1);
  space.assign_dofs();

  Solution Tprev, // previous time step solution, for the time integration method
           Titer; // solution converging during the Newton's iteration

  WeakForm wf(1);
  if(TIME_DISCR == 1) {
    wf.add_biform(0, 0, callback(J_euler), UNSYM, ANY, 1, &Titer);
    wf.add_liform(0, callback(F_euler), ANY, 2, &Titer, &Tprev);
  }
  else {
    wf.add_biform(0, 0, callback(J_cranic), UNSYM, ANY, 1, &Titer);
    wf.add_liform(0, callback(F_cranic), ANY, 2, &Titer, &Tprev);
  }

  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  char title[100];
  ScalarView view("", 0, 0, 600, 600);
  view.fix_scale_width(80);
  ScalarView view2("", 700, 0, 600, 600);
  view2.fix_scale_width(80);

  // setting the Dirichlet lift to be the initial condition
  Titer.set_dirichlet_lift(&space, &pss);
  Tprev.set_dirichlet_lift(&space, &pss);
  nls.set_ic(&Titer, &Titer, PROJ_TYPE);

  // view initial guess for Newton's method
//   sprintf(title, "Initial guess for the Newton's method");
//   view.set_title(title);
//   view.show(&Titer);
//   printf("Click into the image window and press any key.\n");
//   view.wait_for_keypress();

  Solution sln;
  // time stepping
  for(int n = 1; n <= NSTEP; n++)
  {

    info("\n---- Time step %d -----------------------------------------------", n);

    // set initial condition for the Newton's iteration
    // actually needed only when space changes
    // otherwise initial solution vector is that one
    // from the previous time level
    //nls.set_ic(&Titer, &Titer);

    int it = 1;
    double res_l2_norm;
    do
    {
      info("\n---- Time step %d, Newton iter %d ---------------------------------\n", n, it++);

      nls.assemble();
      nls.solve(1, &sln);

      res_l2_norm = nls.get_residuum_l2_norm();
      info("Residuum L2 norm: %g\n", res_l2_norm);
      // want to see Newtons iterations
//       sprintf(title, "Time level %d, Newton iteration %d", n, it-1);
//       view.set_title(title);
//       view.show(&sln);
//       view.wait_for_keypress();

      Titer = sln;

    }
    while (res_l2_norm > NEWTON_TOL);

    // visualization of solution on the n-th time level
    sprintf(title, "Time level %d", n);
    //view.set_min_max_range(90,100);
    view.set_title(title);
    view.show(&Titer);
    //view.wait_for_keypress();

    // uncomment one of the following lines to generate a series of video frames
    //vview.save_numbered_screenshot("velocity%03d.bmp", i, true);
    //pview.save_numbered_screenshot("pressure%03d.bmp", i, true);
    // the frames can then be converted to a video file with the command
    // mencoder "mf://velocity*.bmp" -mf fps=20 -o velocity.avi -ovc lavc -lavcopts vcodec=mpeg4



    // copying result of the Newton's iteration into Tprev
    Tprev.copy(&Titer);
  }

  // wait for keyboard or mouse input
  View::wait("Waiting for all views to be closed.");
  return 0;
}
