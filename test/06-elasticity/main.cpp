#include "hermes2d.h"
#include "solver_umfpack.h"


int bc_types(int marker)
  { return (marker == 8) ? BC_ESSENTIAL : BC_NATURAL; }


const double E = 2e11;
const double nu = 0.3;
const double g = 9.81;
const double rho = 1000;
  
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));
const double l2m = lambda + 2*mu;


scalar bilinear_form_0_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdx_b_dudy_dvdy(l2m, fu, mu, fv, ru, rv); }

scalar bilinear_form_0_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdy_b_dudy_dvdx(lambda, fv, mu, fu, rv, ru); }

scalar bilinear_form_1_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdy_b_dudy_dvdx(lambda, fu, mu, fv, ru, rv); }

scalar bilinear_form_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdx_b_dudy_dvdy(mu, fu, l2m, fv, ru, rv); }
  
scalar linear_form_1(RealFunction* fv, RefMap* rv)
  { return -g * rho * int_v(fv, rv); }


int main(int argc, char* argv[])
{
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  Mesh mesh;
  mesh.load("beam.mesh");
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_towards_vertex(0, 5);
  mesh.refine_towards_vertex(3, 5);
  mesh.refine_all_elements();
  
  H1Space xdisp(&mesh, &shapeset);
  xdisp.set_bc_types(bc_types);
  xdisp.set_uniform_order(6);
  int ndofs = xdisp.assign_dofs();
  
  H1Space ydisp(&mesh, &shapeset);
  ydisp.set_bc_types(bc_types);
  ydisp.set_uniform_order(6);
  ndofs += ydisp.assign_dofs(ndofs);

  WeakForm wf(2);
  wf.add_biform(0, 0, bilinear_form_0_0);
  wf.add_biform(0, 1, bilinear_form_0_1, SYM);
  //wf.add_biform(1, 0, bilinear_form_1_0);
  wf.add_biform(1, 1, bilinear_form_1_1);
  wf.add_liform(1, linear_form_1);
  
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(2, &xdisp, &ydisp);
  sys.set_pss(1, &pss);
  
  Solution xsln, ysln;
  sys.assemble();
  sys.solve(2, &xsln, &ysln);
  
  ScalarView sview("Stress", 0, 250, 1600, 500);
  VonMisesFilter mises(&xsln, &ysln, mu, lambda);
  sview.show(&mises, EPS_NORMAL, FN_VAL_0, &xsln, &ysln, 10000);
  sview.set_min_max_range(0, 1e6);
  
  View::wait();
  return 0;
}

