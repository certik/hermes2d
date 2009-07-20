#include "hermes2d.h"
#include "solver_umfpack.h"

// This example demostrates the adaptive multimesh hp-FEM. The problem
// considered is rooted in linear thermoelasticity: a massive hollow conductor
// is heated by induction and cooled by water running inside. We will approximate the
// x-displacement, y-displacement, and the temperature on individual meshes
// equipped with mutually independent adaptivity mechanisms.
//
// PDE: Linear thermoelasticity
//
// BC: u_1 = u_2 = 0 on Gamma_1
//     du_1/dn = du_2/dn = 0 elsewhere
//     temp = TEMP_INNER on Gamma_4
//     negative heat flux with HEAT_FLUX_OUTER elsewhere
//
// The parameters below can be played with. In particular, compare hp- and
// h-adaptivity via the H_ONLY option, and compare the multi-mesh vs. single-mesh
// using the MULTI parameter.
//

const bool MULTI = true;               // true = use multi-mesh, false = use single-mesh
                                       // Note: in the single mesh option, the meshes are
                                       // forced to be geometrically the same but the
                                       // polynomial degrees can still vary
const bool SAME_ORDERS = false;        // true = when single mesh is used it forces same pol. orders for components
                                       // when multi mesh used, parameter is ignored
const int MESH_REGULARITY = -1;        // specifies level of hanging nodes
                                       // -1 is arbitrary level hangning nodes
                                       // 1, 2, 3,... is 1-irregular mesh, 2-irregular mesh,...
                                       // total regularization (0) is not supported in adaptivity
const int P_INIT_TEMP = 2;             // initial polynomial degrees in temperature mesh
const int P_INIT_DISP = 2;             // initial polynomial degrees for displacement meshes
const double ERR_STOP = 0.02;          // stopping criterion for hp-adaptivity
                                       // (rel. error tolerance between the reference
                                       // and coarse solution in percent)
const double THRESHOLD = 0.3;          // error threshold for element refinement (multi-mesh)
const int STRATEGY = 0;                // refinement strategy (0, 1, 2, 3 - see adapt_h1.cpp for explanation)
const bool H_ONLY = false;             // if H_ONLY == false then full hp-adaptivity takes place, otherwise
                                       // h-adaptivity is used. Use this parameter to check that indeed adaptive
                                       // hp-FEM converges much faster than adaptive h-FEM
const bool ISO_ONLY = false;           // when ISO_ONLY = true, only isometric refinements are done,
                                       // otherwise anisotropic refinements can be taken into account
const int MAX_ORDER = 6;               // maximal order used during adaptivity
const int NDOF_STOP = 40000;           // adaptivity process stops when the number of degrees of freedom grows over
                                       // this limit. This is mainly to prevent h-adaptivity to go on forever.

double HEAT_SRC = 10000.0;       // heat source in the material (caused by induction heating)
double TEMP_INNER = 50;
double HEAT_FLUX_OUTER = -50;
const double E = 2e11;           // steel: E=200 GPa
const double nu = 0.3;
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));
const double l2m = lambda + 2*mu;
const double rho = 8000;
const double g = 9.81;
const double alpha = 13e-6;      // see http://hyperphysics.phy-astr.gsu.edu/hbase/tables/thexp.html

//  Boundary markers:
//  1 - bottom
//  3 - top
//  2 - left & right
//  4 - holes

int bc_types_x(int marker)
  { return (marker == 1) ? BC_ESSENTIAL : BC_NATURAL; }

int bc_types_y(int marker)
  { return (marker == 1) ? BC_ESSENTIAL : BC_NATURAL; }

int bc_types_t(int marker)
  { return (marker == 4) ? BC_ESSENTIAL : BC_NATURAL; }

double bc_values_t(int marker, double x, double y)
  { return (marker == 4) ? TEMP_INNER : HEAT_FLUX_OUTER; }

scalar bilinear_form_unsym_0_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdx_b_dudy_dvdy(l2m, fu, mu, fv, ru, rv); }

scalar bilinear_form_unsym_0_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdy_b_dudy_dvdx(lambda, fv, mu, fu, rv, ru); }

scalar bilinear_form_unsym_0_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return -(3*lambda + 2*mu) * alpha * int_dudx_v(fu, fv, ru, rv); }

scalar bilinear_form_unsym_1_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdy_b_dudy_dvdx(lambda, fu, mu, fv, ru, rv); }

scalar bilinear_form_unsym_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_a_dudx_dvdx_b_dudy_dvdy(mu, fu, l2m, fv, ru, rv); }

scalar bilinear_form_unsym_1_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return -(3*lambda + 2*mu) * alpha * int_dudy_v(fu, fv, ru, rv); }

scalar bilinear_form_unsym_2_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv); }

scalar linear_form_1(RealFunction* fv, RefMap* rv)
  { return -g * rho * int_v(fv, rv); }

scalar linear_form_2(RealFunction* fv, RefMap* rv)
  { return HEAT_SRC * int_v(fv, rv); }

scalar linear_form_surf_2(RealFunction* fv, RefMap* rv, EdgePos* ep)
  { return surf_int_G_v(fv, rv, ep); }


/*
int crit(Element* e)
{
  if (e->is_triangle()) return -1;

  double2 mm[4];
  for (int i = 0; i < 4; i++)
  {
    int j = e->next_vert(i);
    mm[i][0] = (e->vn[i]->x + e->vn[j]->x) / 2;
    mm[i][1] = (e->vn[i]->y + e->vn[j]->y) / 2;
  }

  double l1 = hypot(mm[0][0] - mm[2][0], mm[0][1] - mm[2][1]);
  double l2 = hypot(mm[1][0] - mm[3][0], mm[1][1] - mm[3][1]);

  const double a = 1.2;
  if (l1 > a*l2) return 1;
  if (l2 > a*l1) return 2;

  return -1;
}
*/

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh xmesh, ymesh, tmesh;
  xmesh.load("domain_round_3.mesh"); // master mesh
  ymesh.copy(&xmesh);                // ydisp will share master mesh with xdisp
  tmesh.copy(&xmesh);                // temp will share master mesh with xdisp

  // initialize the shapeset and the cache
  H1ShapesetOrtho shapeset;
  PrecalcShapeset xpss(&shapeset);
  PrecalcShapeset ypss(&shapeset);
  PrecalcShapeset tpss(&shapeset);

  // create the x displacement space
  H1Space xdisp(&xmesh, &shapeset);
  xdisp.set_bc_types(bc_types_x);
  xdisp.set_uniform_order(P_INIT_DISP);

  // create the y displacement space
  H1Space ydisp(MULTI ? &ymesh : &xmesh, &shapeset);
  ydisp.set_bc_types(bc_types_y);
  ydisp.set_uniform_order(P_INIT_DISP);

  // create the temperature space
  H1Space temp(MULTI ? &tmesh : &xmesh, &shapeset);
  temp.set_bc_types(bc_types_t);
  temp.set_bc_values(bc_values_t);
  temp.set_uniform_order(P_INIT_TEMP);

  // initialize the weak formulation
  WeakForm wf(3);
  wf.add_biform(0, 0, bilinear_form_unsym_0_0);
  wf.add_biform(0, 1, bilinear_form_unsym_0_1, SYM);
  wf.add_biform(0, 2, bilinear_form_unsym_0_2);
  wf.add_biform(1, 1, bilinear_form_unsym_1_1);
  wf.add_biform(1, 2, bilinear_form_unsym_1_2);
  wf.add_biform(2, 2, bilinear_form_unsym_2_2);
  wf.add_liform(1, linear_form_1);
  wf.add_liform(2, linear_form_2);
  wf.add_liform_surf(2, linear_form_surf_2);

  // visualization
  OrderView xord("X displacement poly degrees", 0, 0, 850, 400);
  OrderView yord("Y displacement poly degrees", 0, 455, 850, 400);
  OrderView tord("Temperature poly degrees", 0, 885, 850, 400);
  ScalarView sview("Von Mises stress [Pa]", 860, 0, 850, 400);
  ScalarView tview("Temperature [deg C]", 860, 455, 850, 400);
  // disable scales
  //xord.show_scale(false); yord.show_scale(false); tord.show_scale(false);

  GnuplotGraph graph;
  graph.set_captions("", "Degrees of Freedom", "Error (Energy Norm)");
  graph.set_log_y();
  graph.add_row("Refenrence solution", "k", "-", "O");

  GnuplotGraph graph_cpu;
  graph_cpu.set_captions("", "CPU", "error");
  graph_cpu.set_log_y();
  graph_cpu.add_row(MULTI ? "multi-mesh" : "single-mesh", "k", "-", "o");

  Solution xsln, ysln, tsln;
  Solution xrsln, yrsln, trsln;
  UmfpackSolver umfpack;

  int it = 0, ndofs;
  bool done = false;
  double cpu = 0.0;
  do
  {
    printf("\n\n---- it=%d ------------------------------------------------------------------\n\n", it++);
    begin_time();

    ndofs = xdisp.assign_dofs(0);
    ndofs += ydisp.assign_dofs(ndofs);
    ndofs += temp.assign_dofs(ndofs);
    printf("xdof=%d, ydof=%d, tdof=%d\n", xdisp.get_num_dofs(), ydisp.get_num_dofs(), temp.get_num_dofs());

    // solve the coarse problem
    LinSystem ls(&wf, &umfpack);
    ls.set_spaces(3, &xdisp, &ydisp, &temp);
    ls.set_pss(3, &xpss, &ypss, &tpss);
    ls.assemble();
    ls.solve(3, &xsln, &ysln, &tsln);

    cpu += end_time();

    xord.show(&xdisp);
    yord.show(&ydisp);
    tord.show(&temp);

    VonMisesFilter mises(&xsln, &ysln, mu, lambda);
    sview.set_min_max_range(0, 4e9);
    sview.show(&mises, EPS_HIGH);
    tview.show(&tsln, EPS_HIGH);

    // solve the fine (reference) problem
    begin_time();
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(3, &xrsln, &yrsln, &trsln);

    H1OrthoHP hp(3, &xdisp, &ydisp, &temp);
    double error = hp.calc_energy_error_n(3, &xsln, &ysln, &tsln, &xrsln, &yrsln, &trsln,
                   bilinear_form_unsym_0_0, bilinear_form_unsym_0_1, bilinear_form_unsym_0_2,
                   bilinear_form_unsym_1_0, bilinear_form_unsym_1_1, bilinear_form_unsym_1_2,
                   NULL,                    NULL,                    bilinear_form_unsym_2_2) * 100;
    info("\nEstimate of error: %g%%", error);
    graph.add_values(0, xrsln.get_num_dofs() + yrsln.get_num_dofs() + trsln.get_num_dofs(), error);
    graph.save(MULTI ? "conv_m.gp" : "conv_s.gp");
    if (error < ERR_STOP || xdisp.get_num_dofs() + ydisp.get_num_dofs() + temp.get_num_dofs() >= NDOF_STOP) done = true;
    else hp.adapt(THRESHOLD, STRATEGY, H_ONLY, ISO_ONLY, MESH_REGULARITY, MAX_ORDER, SAME_ORDERS);

    cpu += end_time();
    graph_cpu.add_values(0, cpu, error);
    graph_cpu.save(MULTI ? "cpu_m.gp" : "cpu_s.gp");
  }
  while (!done);

  // wait for keypress or mouse input
  View::wait();
  return 0;
};

