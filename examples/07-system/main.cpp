#include "hermes2d.h"
#include "solver_umfpack.h"

const double E  = 200e9; // steel: 200GPa
const double nu = 0.3;
const double f  = 1e4;   // 10^4 N
  
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));


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
  xdisp.set_uniform_order(8);
  int ndofs = xdisp.assign_dofs(0);

  // create the y displacement space
  H1Space ydisp(&mesh, &shapeset);
  ydisp.set_bc_types(bc_types_y);
  ydisp.set_bc_values(bc_values_y);
  ydisp.set_uniform_order(8);
  ndofs += ydisp.assign_dofs(ndofs);

  // initialize the weak formulation
  WeakForm wf(2);
  wf.add_biform(0, 0, bilinear_form_0_0, SYM);
  wf.add_biform(0, 1, bilinear_form_0_1, SYM);
  wf.add_biform(1, 1, bilinear_form_1_1, SYM);
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
  VonMisesFilter stress(&xsln, &ysln, mu, lambda);
  view.show(&stress, EPS_HIGH, FN_VAL_0, &xsln, &ysln, 1.5e5);
  
  View::wait();
  return 0;
}
