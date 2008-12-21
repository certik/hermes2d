#include "hermes2d.h"
#include "solver_umfpack.h"


const double Re = 1000;  // Reynolds number
const double tau = 0.05; // time step


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


scalar bilinear_form_unsym_0_0_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv) / Re +
           int_u_v(fu, fv, ru, rv) / tau +
           int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv); }

scalar bilinear_form_unsym_0_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return -int_u_dvdx(fu, fv, ru, rv); }

scalar bilinear_form_unsym_1_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return -int_u_dvdy(fu, fv, ru, rv); }

scalar bilinear_form_unsym_2_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_dudx_v(fu, fv, ru, rv); }

scalar bilinear_form_unsym_2_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_dudy_v(fu, fv, ru, rv); }

scalar linear_form_0(RealFunction* fv, RefMap* rv)
  { return int_u_v(&xprev, fv, rv, rv) / tau; }

scalar linear_form_1(RealFunction* fv, RefMap* rv)
  { return int_u_v(&yprev, fv, rv, rv) / tau; }


int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);

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

  // zero initial xprev and yprev values
  xprev.set_zero(&mesh);
  yprev.set_zero(&mesh);

  // set up weak formulation
  WeakForm wf(3);
  wf.add_biform(0, 0, bilinear_form_unsym_0_0_1_1, UNSYM, 0, 2, &xprev, &yprev);
  wf.add_biform(0, 2, bilinear_form_unsym_0_2);
  wf.add_biform(1, 1, bilinear_form_unsym_0_0_1_1, UNSYM, 0, 2, &xprev, &yprev);
  wf.add_biform(1, 2, bilinear_form_unsym_1_2);
  wf.add_biform(2, 0, bilinear_form_unsym_2_0);
  wf.add_biform(2, 1, bilinear_form_unsym_2_1);
  wf.add_liform(0, linear_form_0, 0, 1, &xprev);
  wf.add_liform(1, linear_form_1, 0, 1, &yprev);

  // visualization
  VectorView vview("velocity [m/s]", 0, 0, 1200, 470);
  ScalarView pview("pressure [Pa]", 0, 500, 1200, 470);
  vview.set_min_max_range(0, 2.0);
  //vview.show_scale(false);
  //pview.show_scale(false);
  pview.show_mesh(false);
  //pview.set_min_max_range(-1, 1);
  //pview.set_num_palette_steps(2);

  // set up solver and linear system
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(3, &xvel, &yvel, &press);
  sys.set_pss(1, &pss);

  // main loop
  for (int i = 0; i < 1000; i++)
  {
    printf("\n*** Iteration %d ***\n", i);

    // assemble and solve
    Solution xsln, ysln, psln;
    sys.assemble();
    sys.solve(3, &xsln, &ysln, &psln);

    // visualization
    vview.show(&xsln, &ysln, 2*EPS_LOW);
    pview.show(&psln);

    // uncomment one of the following lines to generate a series of video frames
    //vview.save_numbered_screenshot("velocity%03d.bmp", i, true);
    //pview.save_numbered_screenshot("pressure%03d.bmp", i, true);
    // the frames can then be converted to a video file with the command
    // mencoder "mf://velocity*.bmp" -mf fps=20 -o velocity.avi -ovc lavc -lavcopts vcodec=mpeg4

    xprev = xsln;
    yprev = ysln;
  }

  // done
  hermes2d_finalize();
}
