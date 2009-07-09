#include "hermes2d.h"
#include "solver_umfpack.h"

// The time-dependent laminar incompressible Navier-Stokes equations are
// discretized in time via the implicit Euler method. The convective term
// is linearized simply by replacing the velocity in front of the nabla
// operator with the velocity from last time step.
//
// PDE: incompressible Navier-Stokes equations in the form
// \partial v / \partial t - \Delta v / Re + (v \cdot \nabla) v + \nabla p = 0,
// div v = 0
//
// BC: u_1 is a time-dependent constant and u_2 = 0 on Gamma_4 (inlet)
//     u_1 = u_2 = 0 on Gamma_1 (bottom), Gamma_3 (top) and Gamma_5 (obstacle)
//     "do nothing" on Gamma_2 (outlet)
//
// TODO: Implement Crank-Nicolson so that comparisons with implicit Euler can be made
//
// The following parameters can be changed:
//

double RE = 200.0;             // Reynolds number
double VEL_INLET = 1.0;        // inlet velocity (reached after STARTUP_TIME)
double STARTUP_TIME = 1.0;     // during this time, inlet velocity increases gradually
                               // from 0 to VEL_INLET, then it stays constant
double TAU = 0.1;              // time step
double FINAL_TIME = 3000.0;    // length of time interval
int P_INIT_VEL = 2;            // initial polynomial degree for velocity components
int P_INIT_PRESSURE = 1;       // initial polynomial degree for pressure
                               // Note: P_INIT_VEL should always be greater than
                               // P_INIT_PRESSURE because of the inf-sup condition
double H = 5;                  // domain height (necessary to define the parabolic
                               // velocity profile at inlet)

//  to better understand boundary conditions
int marker_bottom = 1;
int marker_right  = 2;
int marker_top = 3;
int marker_left = 4;
int marker_obstacle = 5;

// global time variable
double TIME = 0;

// definition of boundary conditions
int xvel_bc_type(int marker) {
  if (marker == 2) return BC_NONE;
  else return BC_ESSENTIAL;
}

int yvel_bc_type(int marker) {
  if (marker == 2) return BC_NONE;
  else return BC_ESSENTIAL;
}

int press_bc_type(int marker)
  { return BC_NONE; }

scalar xvel_bc_value(int marker, double x, double y) {
  if (marker == 4) {
    // time-dependent inlet velocity
    //double val_y = VEL_INLET; //constant profile
    double val_y = VEL_INLET * y*(H-y) / (H/2.)/(H/2.); //parabolic profile with peak VEL_INLET at y = H/2
    if (TIME <= STARTUP_TIME) return val_y * TIME/STARTUP_TIME;
    else return val_y;
  }
  else return 0;
}


// velocities from the previous time step
Solution xprev, yprev;

scalar bilinear_form_sym_0_0_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv) / RE +
           int_u_v(fu, fv, ru, rv) / TAU; }

scalar bilinear_form_unsym_0_0_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv); }

scalar bilinear_form_unsym_0_2(RealFunction* fp, RealFunction* fv, RefMap* rp, RefMap* rv)
  { return -int_u_dvdx(fp, fv, rp, rv); }

scalar bilinear_form_unsym_1_2(RealFunction* fp, RealFunction* fv, RefMap* rp, RefMap* rv)
  { return -int_u_dvdy(fp, fv, rp, rv); }

scalar linear_form_0(RealFunction* fv, RefMap* rv)
  { return int_u_v(&xprev, fv, xprev.get_refmap(), rv) / TAU; }

scalar linear_form_1(RealFunction* fv, RefMap* rv)
  { return int_u_v(&yprev, fv, yprev.get_refmap(), rv) / TAU; }


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");

  // a-priori mesh refinements
  mesh.refine_all_elements();
  mesh.refine_towards_boundary(5, 4, false);
  mesh.refine_towards_boundary(1, 4);
  mesh.refine_towards_boundary(3, 4);

  // display the mesh
  //MeshView mview("Hello world!", 100, 100, 1100, 400);
  //mview.show(&mesh);
  //mview.wait_for_keypress();

  // initialize the shapesets and the cache
  H1ShapesetBeuchler shapeset_v;
  PrecalcShapeset pss_v(&shapeset_v);
  L2Shapeset shapeset_p;
  PrecalcShapeset pss_p(&shapeset_p);

  // H1 spaces for velocities and L2 for pressure
  H1Space xvel(&mesh, &shapeset_v);
  H1Space yvel(&mesh, &shapeset_v);
  //H1Space press(&mesh, &shapeset);
  L2Space press(&mesh, &shapeset_p);

  // initialize boundary conditions
  xvel.set_bc_types(xvel_bc_type);
  xvel.set_bc_values(xvel_bc_value);
  yvel.set_bc_types(yvel_bc_type);
  press.set_bc_types(press_bc_type);

  // set velocity and pressure polynomial degrees
  xvel.set_uniform_order(P_INIT_VEL);
  yvel.set_uniform_order(P_INIT_VEL);
  press.set_uniform_order(P_INIT_PRESSURE);

  // assign degrees of freedom
  int ndofs = 0;
  ndofs += xvel.assign_dofs(ndofs);
  ndofs += yvel.assign_dofs(ndofs);
  ndofs += press.assign_dofs(ndofs);

  // initial BC: xprev and yprev are zero
  xprev.set_zero(&mesh);
  yprev.set_zero(&mesh);

  // set up weak formulation
  WeakForm wf(3);
  wf.add_biform(0, 0, bilinear_form_sym_0_0_1_1, SYM);
  wf.add_biform(0, 0, bilinear_form_unsym_0_0_1_1, UNSYM, ANY, 2, &xprev, &yprev);
  wf.add_biform(1, 1, bilinear_form_sym_0_0_1_1, SYM);
  wf.add_biform(1, 1, bilinear_form_unsym_0_0_1_1, UNSYM, ANY, 2, &xprev, &yprev);
  wf.add_biform(0, 2, bilinear_form_unsym_0_2, ANTISYM);
  wf.add_biform(1, 2, bilinear_form_unsym_1_2, ANTISYM);
  wf.add_liform(0, linear_form_0, ANY, 1, &xprev);
  wf.add_liform(1, linear_form_1, ANY, 1, &yprev);

  // visualization
  VectorView vview("velocity [m/s]", 0, 0, 1500, 470);
  ScalarView pview("pressure [Pa]", 0, 530, 1500, 470);
  vview.set_min_max_range(0, 1.6);
  //pview.set_min_max_range(-0.9, 0.9);
  pview.show_mesh(false);

  // set up the linear system
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(3, &xvel, &yvel, &press);
  sys.set_pss(3, &pss_v, &pss_v, &pss_p);

  // main loop
  char title[100];
  int num_time_steps = FINAL_TIME / TAU;
  for (int i = 1; i <= num_time_steps; i++)
  {
    TIME += TAU;

    info("\n---- Time step %d, time = %g -----------------------------------", i, TIME);

    // this is needed to update the time-dependent boundary conditions
    ndofs = 0;
    ndofs += xvel.assign_dofs(ndofs);
    ndofs += yvel.assign_dofs(ndofs);
    ndofs += press.assign_dofs(ndofs);

    // assemble and solve
    Solution xsln, ysln, psln;
    psln.set_zero(&mesh);
    sys.assemble();
    sys.solve(3, &xsln, &ysln, &psln);

    // visualization
    sprintf(title, "Velocity, time %g", TIME);
    vview.set_title(title);
    vview.show(&xprev, &yprev, EPS_LOW);
    sprintf(title, "Pressure, time %g", TIME);
    pview.set_title(title);
    pview.show(&psln);

    // uncomment one of the following lines to generate a series of video frames
    //vview.save_numbered_screenshot("velocity%03d.bmp", i, true);
    //pview.save_numbered_screenshot("pressure%03d.bmp", i, true);
    // the frames can then be converted to a video file with the command
    // mencoder "mf://velocity*.bmp" -mf fps=20 -o velocity.avi -ovc lavc -lavcopts vcodec=mpeg4

    xprev = xsln;
    yprev = ysln;
  }

  View::wait();
}
