#include "hermes2d.h"
#include "solver_umfpack.h"

const double Re = 700;  // Reynolds number  FIXME - char. length is not 1
const double tau = 0.05; // time step

const int marker_bottom = 1;
const int marker_right  = 2;


// definition of boundary conditions
int xvel_bc_type(int marker)
  { return (marker != 2) ? BC_ESSENTIAL : BC_NONE; }

scalar xvel_bc_value(int marker, double x, double y)
  { return (marker != 5) ? 1 : 0; }

int yvel_bc_type(int marker)
  { return (marker != 2) ? BC_ESSENTIAL : BC_NONE; }

int press_bc_type(int marker)
  { return BC_NONE; }


// velocities from the previous time step
Solution xprev, yprev;

scalar bilinear_form_sym_0_0_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv) / Re +
           int_u_v(fu, fv, ru, rv) / tau; }

scalar bilinear_form_unsym_0_0_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv); }

scalar bilinear_form_unsym_0_2(RealFunction* fp, RealFunction* fv, RefMap* rp, RefMap* rv)
  { return -int_u_dvdx(fp, fv, rp, rv); }

scalar bilinear_form_unsym_1_2(RealFunction* fp, RealFunction* fv, RefMap* rp, RefMap* rv)
  { return -int_u_dvdy(fp, fv, rp, rv); }

scalar linear_form_0(RealFunction* fv, RefMap* rv)
  { return int_u_v(&xprev, fv, xprev.get_refmap(), rv) / tau; }

scalar linear_form_1(RealFunction* fv, RefMap* rv)
  { return int_u_v(&yprev, fv, yprev.get_refmap(), rv) / tau; }


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("cylinder4.mesh");
  mesh.refine_towards_boundary(5, 3);

  H1ShapesetBeuchler shapeset;
  PrecalcShapeset pss(&shapeset);

  // spaces for velocities and pressure
  H1Space xvel(&mesh, &shapeset);
  H1Space yvel(&mesh, &shapeset);
  H1Space press(&mesh, &shapeset);

  // initialize boundary conditions
  xvel.set_bc_types(xvel_bc_type);
  xvel.set_bc_values(xvel_bc_value);
  yvel.set_bc_types(yvel_bc_type);
  press.set_bc_types(press_bc_type);

  // set velocity and pressure polynomial degrees
  xvel.set_uniform_order(2);
  yvel.set_uniform_order(2);
  press.set_uniform_order(1);

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
  VectorView vview("velocity [m/s]", 0, 0, 1200, 470);
  ScalarView pview("pressure [Pa]", 0, 500, 1200, 470);
  vview.set_min_max_range(0, 1.9);
  pview.set_min_max_range(-0.9, 0.9);
  pview.show_mesh(false);

  // set up the linear system
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(3, &xvel, &yvel, &press);
  sys.set_pss(1, &pss);

  // main loop
  for (int i = 0; i < 1000; i++)
  {
    info("\n*** Iteration %d ***", i);
    
    // assemble and solve
    Solution xsln, ysln, psln;
    sys.assemble();
    sys.solve(3, &xsln, &ysln, &psln);

    // visualization
    vview.show(&xprev, &yprev, EPS_LOW);
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
