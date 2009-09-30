#include "hermes2d.h"
#include "solver_umfpack.h"

//
//  This example demonstrates the employment of the Newton's method for
//  a nonlinear complex-valued time-dependent PDE (the Gross-Pitaevski
//  equation describing the behavior of Einstein-Bose quantum gases)
//  discretized implicitly in time (via implicit Euler or Crank-Nicolson).
//  Some problem parameters can be changed below.
//
//  PDE: non-stationary complex Gross-Pitaevski equation
//  describing resonances in Bose-Einstein condensates
//
//  ih \partial \psi/\partial t = -h^2/(2m) \Delta \psi +
//  g \psi |\psi|^2 + 1/2 m \omega^2 (x^2 + y^2) \psi
//
//  Domain: square of edge length 2 centered at origin
//
//  BC:  homogeneous Dirichlet everywhere on th eboundary
//
//  Time-stepping: either implicit Euler or Crank-Nicolson
//

/********** Problem parameters ***********/

double H = 1;              //Planck constant 6.626068e-34;
double M = 1;
double G = 1;
double OMEGA = 1;
int TIME_DISCR = 2;        // 1 for implicit Euler, 2 for Crank-Nicolson
int PROJ_TYPE = 2;         // 1 for H1 projections, 2 for L2 projections
double T_FINAL = 2;        // time interval length
double TAU = 0.001;        // time step
int P_INIT = 2;            // initial polynomial degree
int REF_INIT = 3;          // number of initial uniform refinements
double NEWTON_TOL = 1e-3;  // convergence criterion for the Newton's method

/********** Definition of initial conditions ***********/

scalar fn_init(double x, double y, scalar& dx, scalar& dy) {
  return complex(exp(-10*(x*x + y*y)), 0);
}

/********** Definition of boundary conditions ***********/

int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

complex bc_values(int marker, double x, double y)
{
  return complex(0.0, 0.0);
}

/********** Definition of Jacobian matrices and residual vectors ***********/

# include "forms.cpp"

/********** Definition of linear and bilinear forms for Hermes ***********/

Solution Psi_prev, // previous time step solution, for the time integration method
         Psi_iter; // solution converging during the Newton's iteration

// Implicit Euler method (1st-order in time)
complex linear_form_0_euler(RealFunction* fv, RefMap* rv)
{ return F_euler(&Psi_prev, &Psi_iter, fv, rv); }
complex bilinear_form_0_0_euler(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ return J_euler(&Psi_iter, fu, fv, ru, rv); }

// Crank-Nicolson method (2nd-order in time)
complex linear_form_0_cranic(RealFunction* fv, RefMap* rv)
{ return F_cranic(&Psi_prev, &Psi_iter, fv, rv); }
complex bilinear_form_0_0_cranic(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ return J_cranic(&Psi_iter, fu, fv, ru, rv); }


// *************************************************************

int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("square.mesh");
  for(int i = 0; i < REF_INIT; i++) mesh.refine_all_elements();

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
  space.assign_dofs();

  WeakForm wf(1);
  if(TIME_DISCR == 1) {
    wf.add_biform(0, 0, bilinear_form_0_0_euler, UNSYM, ANY, 1, &Psi_iter);
    wf.add_liform(0, linear_form_0_euler, ANY, 2, &Psi_iter, &Psi_prev);
  }
  else {
    wf.add_biform(0, 0, bilinear_form_0_0_cranic, UNSYM, ANY, 1, &Psi_iter);
    wf.add_liform(0, linear_form_0_cranic, ANY, 2, &Psi_iter, &Psi_prev);
  }

  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  char title[100];
  ScalarView view("", 0, 0, 700, 600);
  //view.set_min_max_range(-0.5,0.5);

  // setting initial condition at zero time level
  Psi_prev.set_exact(&mesh, fn_init);
  nls.set_ic(&Psi_prev, &Psi_prev, PROJ_TYPE);
  Psi_iter.copy(&Psi_prev);

  // view initial guess for Newton's method
  /*
  sprintf(title, "Initial guess for the Newton's method");
  view.set_title(title);
  view.show(&Psi_iter);
  view.wait_for_keypress();
  */

  Solution sln;
  // time stepping
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int n = 1; n <= nstep; n++)
  {

    info("\n---- Time step %d -----------------------------------------------", n);

    // set initial condition for the Newton's iteration
    // actually needed only when space changes
    // otherwise initial solution vector is that one
    // from the previous time level
    //nls.set_ic(&Psi_iter, &Psi_iter);

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
      /*
       sprintf(title, "Time level %d, Newton iteration %d", n, it-1);
       view.set_title(title);
       view.show(&sln);
       view.wait_for_keypress();
      */

      Psi_iter = sln;

    }
    while (res_l2_norm > NEWTON_TOL);

    // visualization of solution on the n-th time level
    sprintf(title, "Time level %d", n);
    //view.set_min_max_range(90,100);
    view.set_title(title);
    view.show(&Psi_iter);
    //view.wait_for_keypress();

    // uncomment one of the following lines to generate a series of video frames
    //view.save_numbered_screenshot("sol%03d.bmp", n, true);
    //pview.save_numbered_screenshot("pressure%03d.bmp", i, true);
    // the frames can then be converted to a video file with the command
    // mencoder "mf://velocity*.bmp" -mf fps=20 -o velocity.avi -ovc lavc -lavcopts vcodec=mpeg4



    // copying result of the Newton's iteration into Psi_prev
    Psi_prev.copy(&Psi_iter);
  }

  printf("Click into the image window and press 'q' to finish.\n");
  View::wait();
  return 0;
}
