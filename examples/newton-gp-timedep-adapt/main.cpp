#include "hermes2d.h"
#include "solver_umfpack.h"

//
//  This example shows how to combine automatic adaptivity with the Newton's
//  method for a nonlinear complex-valued time-dependent PDE (the Gross-Pitaevski
//  equation describing the behavior of Einstein-Bose quantum gases)
//  discretized implicitly in time (via implicit Euler or Crank-Nicolson).
//  Some problem parameters can be changed below.
//
//  PDE: non-stationary complex Gross-Pitaevski equation
//  describing resonances in Bose-Einstein condensates
//
//  hp-adaptivity with dynamical meshes
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

double H = 1;         //Planck constant 6.626068e-34;
double M = 1;
double G = 1;
double OMEGA = 1;
int TIME_DISCR = 2;                // 1 for implicit Euler, 2 for Crank-Nicolson
int PROJ_TYPE = 1;                 // 1 for H1 projections, 0 for L2 projections
double T_FINAL = 200.0;            // time interval length
double TAU = 0.005;                // time step
int P_INIT = 1;                    // initial polynomial degree
int REF_INIT = 2;                  // number of initial uniform refinements

// Newton parameters
double NEWTON_TOL_COARSE = 0.05;   // stopping criterion for Newton on coarse mesh
double NEWTON_TOL_REF = 0.5;       // stopping criterion for Newton on fine mesh
                                   // (the ref. solution does not to be super accurate)

// adaptivity parameters
int UNREF_FREQ = 1;                // every UNREF_FREQth time step the mesh
                                   // is unrefined
double SPACE_L2_TOL = 1.0;         // stopping criterion for hp-adaptivity
                                   // (relative error between reference and coarse solution in percent)
double THR = 0.3;                  // parameter of the adaptivity procedure indicating how many
                                   // refinements will be done per adaptivity step
int STRATEGY = 1;                  // adaptive strategy (0, 1, 2 or 3)
bool H_ONLY = false;               // only perform h-adaptivity
bool ISO_ONLY = false;             // only consider isotropic refinements
int MAX_P = 5;                     // maximum polynomial order allowed in hp-adaptivity
                                   // had to be limited due to complicated integrals
int SHOW_MESHES_IN_TIME_STEP = 0;  // if nonzero, all meshes during every time step are
                                   // shown, else only meshes at the end of every time step.

/********** Definition of initial conditions ***********/

scalar fn_init(double x, double y, scalar& dx, scalar& dy) {
  return exp(-20*(x*x + y*y));
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
  Mesh mesh, basemesh;
  basemesh.load("square.mesh");
  for(int i = 0; i < REF_INIT; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
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
  ScalarView view("", 0, 0, 600, 600);
  //ScalarView realview("", 0, 0, 600, 600);
  //ScalarView imagview("", 700, 0, 600, 600);
  ScalarView magview("", 0, 0, 600, 600);
  //ScalarView angleview("", 700, 700, 600, 600);
  OrderView ordview("", 700, 0, 600, 600);

  GnuplotGraph graph_err;
  graph_err.set_captions("","Time step","Error");
  graph_err.add_row();
  GnuplotGraph graph_dofs;
  graph_dofs.set_captions("","Time step","DOFs");
  graph_dofs.add_row();

  // setting initial condition at zero time level
  Psi_prev.set_exact(&mesh, fn_init);
  Psi_iter.set_exact(&mesh, fn_init);
  nls.set_ic(&Psi_iter, &Psi_iter, PROJ_TYPE);

  // showing initial condition
  //sprintf(title, "Initial condititon");
  //view.set_title(title);
  //view.show(&Psi_prev);
  //ordview.show(&space);
  //view.wait_for_keypress(); // this may cause graphics problems

  Solution sln, rsln;
  // time stepping
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int n = 1; n <= nstep; n++)
  {

    info("\n---- Time step %d -----------------------------------------------------------------", n);

    // unrefinements
    if (n % UNREF_FREQ == 0) {              // frequency of unrefinements
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);
    }

    int at = 0;
    bool done = false;
    double err;
    do
    {
     info("\n---- Time step %d, adaptivity step %d ---------------------------------------------\n", n, ++at);


      int it = 1;
      double res_l2_norm;
      space.assign_dofs();
      if (n > 1 || at > 1) nls.set_ic(&rsln, &Psi_iter);
      else nls.set_ic(&Psi_iter, &Psi_iter);
      do
      {
        info("\n---- Time step %d, adaptivity step %d, Newton step %d (coarse solution) ----------\n", n, at, it++);

        nls.assemble();
        nls.solve(1, &sln);

        res_l2_norm = nls.get_residuum_l2_norm();
        info("Residuum L2 norm: %g\n", res_l2_norm);

        Psi_iter.copy(&sln);
      }
      while (res_l2_norm > NEWTON_TOL_COARSE);

      // reference solution
      it = 1;
      RefNonlinSystem rs(&nls);
      rs.prepare();
      if (n > 1 || at > 1) rs.set_ic(&rsln, &Psi_iter);
      else rs.set_ic(&Psi_iter, &Psi_iter);
      do
      {
        info("\n---- Time step %d, adaptivity step %d, Newton step %d (reference solution) --------\n", n, at, it++);

        rs.assemble();
        rs.solve(1, &rsln);

        res_l2_norm = rs.get_residuum_l2_norm();
        info("Residuum L2 norm: %g\n", res_l2_norm);

        Psi_iter.copy(&rsln);
      }
      while (res_l2_norm > NEWTON_TOL_REF);

      // visualization of intermediate solution
      // and mesh during adaptivity
      if(SHOW_MESHES_IN_TIME_STEP) {
        sprintf(title, "Magnitude, time level %d", n);
        magview.set_title(title);
        AbsFilter mag(&sln);
        magview.show(&mag);            // to see the magnitude of the solution
        sprintf(title, "hp-mesh, time level %d", n);
        ordview.set_title(title);
        ordview.show(&space);          // to see hp-mesh during the process
      }

      H1OrthoHP hp(1, &space);
      err = hp.calc_error(&sln, &rsln) * 100;   // relative l2-error in percent
      info("Error: %g", err);
      if (err < SPACE_L2_TOL) done = true;
      else hp.adapt(THR, STRATEGY, H_ONLY, ISO_ONLY, MAX_P);

    }
    while (!done);

    // showing real part of the solution
    //sprintf(title, "Real Component, time level %d", n);
    //realview.set_title(title);
    //RealFilter real(&sln);
    //realview.show(&real);       // to see the solution - real component

    // showing imaginary part of the solution
    //sprintf(title, "Imag Component, time level %d", n);
    //imagview.set_title(title);
    //ImagFilter imag(&sln);
    //imagview.show(&imag);      // to see the solution - imag component

    // showing magnitude of the solution
    sprintf(title, "Magnitude, time level %d", n);
    magview.set_title(title);
    AbsFilter mag(&sln);
    magview.show(&mag);      // to see the magnitude of the solution

    // showing the angle of the solution
    //sprintf(title, "Angle, Time level %d", n);
    //angleview.set_title(title);
    //AngleFilter angle(&sln);
    //angleview.show(&angle);      // to see the angle of the solution

    // showing hp-mesh
    sprintf(title, "hp-mesh, time level %d", n);
    ordview.set_title(title);
    ordview.show(&space);      // to see hp-mesh

    // to save solutions
    //ordview.save_numbered_screenshot("./video/mesh%04d.bmp", n, true);
    //magview.save_numbered_screenshot("./video/sol%04d.bmp", n, true);

    graph_err.add_values(0, n, err);
    graph_err.save("error.txt");
    graph_dofs.add_values(0, n, space.get_num_dofs());
    graph_dofs.save("dofs.txt");

    // copying result of the Newton's iteration into Psi_prev
    Psi_prev.copy(&rsln);
    printf("hello\n");
  }

  printf("Click into the image window and press 'q' to finish.\n");
  View::wait();
  return 0;
}
