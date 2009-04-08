#include "hermes2d.h"
#include "solver_umfpack.h"

const double tol = 0.01; // error tolerance in percent
const double thr = 0.3;  // error threshold for element refinement

const double E  = 200e9; // Young modulus for steel: 200GPa
const double nu = 0.3;   // Poisson ratio
const double f  = 1e3;   // load force: 10^3 N
  
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));

const int marker_left = 1;
const int marker_top = 2;


int bc_types_xy(int marker)
  { return (marker == marker_left) ? BC_ESSENTIAL : BC_NATURAL; }


scalar bilinear_form_0_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdx_b_dudy_dvdy(lambda+2*mu, fu, mu, fv, ru, rv); }

scalar bilinear_form_0_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdy_b_dudy_dvdx(lambda, fv, mu, fu, rv, ru); }

scalar bilinear_form_1_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdy_b_dudy_dvdx(lambda, fu, mu, fv, ru, rv); }

/*scalar bilinear_form_0_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return (lambda + mu) * int_dudy_dvdx(fu, fv, ru, rv); }

scalar bilinear_form_1_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return (lambda + mu) * int_dudx_dvdy(fu, fv, ru, rv); }*/
  
scalar bilinear_form_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdx_b_dudy_dvdy(mu, fu, lambda+2*mu, fv, ru, rv); }

scalar linear_form_1_surf_top(RealFunction* fv, RefMap* rv, EdgePos* ep)
  { return -f * surf_int_v(fv, rv, ep); }


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh xmesh, ymesh;
  xmesh.load("bracket.mesh");
  ymesh.copy(&xmesh);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset xpss(&shapeset);
  PrecalcShapeset ypss(&shapeset); // fixme: this shouldn't be necessary

  // create the x displacement space
  H1Space xdisp(&xmesh, &shapeset);
  xdisp.set_bc_types(bc_types_xy);
  xdisp.set_uniform_order(1);

  // create the y displacement space
  H1Space ydisp(&ymesh, &shapeset);
  ydisp.set_bc_types(bc_types_xy);
  ydisp.set_uniform_order(1);

  // initialize the weak formulation
  WeakForm wf(2);
  wf.add_biform(0, 0, bilinear_form_0_0, SYM);
  wf.add_biform(0, 1, bilinear_form_0_1, SYM);
  wf.add_biform(1, 1, bilinear_form_1_1, SYM);
  wf.add_liform_surf(1, linear_form_1_surf_top, marker_top);
  
  ScalarView sview("Von Mises stress [Pa]", 0, 300, 800, 800);
  OrderView  xoview("X polynomial orders", 0, 0, 800, 800);
  OrderView  yoview("Y polynomial orders", 810, 0, 800, 800);

  GnuplotGraph graph;
  graph.set_captions("Error Convergence", "Degrees of Freedom", "Error [%]");
  graph.add_row("hp-adaptivity", "k", "-", "o");
  graph.set_log_y();

  Solution xsln, ysln;
  Solution xrsln, yrsln;
  UmfpackSolver umfpack;

  char filename[200];
  
  int it = 1;
  while (1)
  {
    info("\n---- Iteration %d ---------------------------------------------\n", it++);
	  
    int ndofs = xdisp.assign_dofs();
    ndofs += ydisp.assign_dofs(ndofs);
	  
    // solve the coarse problem
    LinSystem ls(&wf, &umfpack);
    ls.set_spaces(2, &xdisp, &ydisp);
    ls.set_pss(2, &xpss, &ypss);
    ls.assemble();
    ls.solve(2, &xsln, &ysln);
  
    // visualize the solution
    VonMisesFilter stress(&xsln, &ysln, mu, lambda);
    sview.set_min_max_range(0, 3e4);
    sview.show(&stress);
    xoview.show(&xdisp);
    yoview.show(&ydisp);
    info("");
  
    // solve the fine (reference) problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(2, &xrsln, &yrsln);
    
    // calculate errors and adapt the solution
    H1OrthoHP hp(2, &xdisp, &ydisp);
    double error = hp.calc_energy_error_2(&xsln, &ysln, &xrsln, &yrsln, 
                                          bilinear_form_0_0, bilinear_form_0_1,
                                          bilinear_form_1_0, bilinear_form_1_1) * 100;
    if (error < tol) break;
    hp.adapt(thr);

    graph.add_values(0, xdisp.get_num_dofs() + ydisp.get_num_dofs(), error);
    graph.save("conv.txt");
    
  }
  
  View::wait();
  return 0;
}
