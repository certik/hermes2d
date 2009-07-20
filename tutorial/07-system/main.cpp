#include "hermes2d.h"
#include "solver_umfpack.h"

// This example explains how to create two spaces over a mesh and use them
// to solve a simple problem of linear elasticity. At the end, VonMises
// filter is used to visualize the stress.
//
// PDE: Lame equations of linear elasticity
//
// BC: du_1/dn = f_1 on Gamma_3 and 0 elsewhere
//     du_2/dn = f_2 on Gamma_3 and 0 Gamma_2, Gamma_4, Gamma_5
//     u_2 = 0 on Gamma_1
//
// The following parameters can be changed:
//

const double E  = 200e9;                                   // Young modulus (steel)
const double nu = 0.3;                                     // Poisson ratio
const double f  = 1e4;                                     // external force

const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));  // Lame constant
const double mu = E / (2*(1 + nu));                        // Lame constant

int P_INIT = 8;

int bc_types_x(int marker)
  { return BC_NATURAL; }

int bc_types_y(int marker)
  { return (marker == 1) ? BC_ESSENTIAL : BC_NATURAL; }

double bc_values_y(EdgePos* ep)
  { return (ep->marker == 3) ? f : 0.0; }


scalar bilinear_form_0_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdx_b_dudy_dvdy(lambda+2*mu, fu, mu, fv, ru, rv); }

scalar bilinear_form_0_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdy_b_dudy_dvdx(lambda, fv, mu, fu, rv, ru); }

scalar bilinear_form_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdx_b_dudy_dvdy(mu, fu, lambda+2*mu, fv, ru, rv); }

scalar linear_form_1_surf(RealFunction* fv, RefMap* rv, EdgePos* ep)
  { return surf_int_G_v(fv, rv, ep); }


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("sample.mesh");

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create the x displacement space
  H1Space xdisp(&mesh, &shapeset);
  xdisp.set_bc_types(bc_types_x);
  xdisp.set_uniform_order(P_INIT);
  int ndofs = xdisp.assign_dofs(0);

  // create the y displacement space
  H1Space ydisp(&mesh, &shapeset);
  ydisp.set_bc_types(bc_types_y);
  ydisp.set_bc_values(bc_values_y);
  ydisp.set_uniform_order(8);
  ndofs += ydisp.assign_dofs(ndofs);

  // initialize the weak formulation
  WeakForm wf(2);
  wf.add_biform(0, 0, bilinear_form_0_0, SYM);  // note that only one symmetric part is
  wf.add_biform(0, 1, bilinear_form_0_1, SYM);  // added in the case of symmetric bilinear
  wf.add_biform(1, 1, bilinear_form_1_1, SYM);  // forms
  wf.add_liform_surf(1, linear_form_1_surf);

  // initialize the linear system and solver
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(2, &xdisp, &ydisp);
  sys.set_pss(1, &pss);

  // assemble the stiffness matrix and solve the system
  Solution xsln, ysln;
  sys.assemble();
  sys.solve(2, &xsln, &ysln);

  // visualize the solution
  ScalarView view("Von Mises stress [Pa]", 50, 50, 1200, 600);
  VonMisesFilter stress(&xsln, &ysln, lambda, mu);
  view.show(&stress, EPS_HIGH, FN_VAL_0, &xsln, &ysln, 1.5e5);

  View::wait();
  return 0;
}
