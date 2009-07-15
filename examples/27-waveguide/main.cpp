#include "hermes2d.h"
#include "solver_umfpack.h"

// This example solved adaptively the electric field in a simplified microwave oven.
// The waves are generated using a harmonic surface current on the right-most edge,
// inside a small square cavity attached to the left edge of the main oven cavity.
// (Such small cavity is present in every microwave oven). There is a circular
// load located in the middle of the main cavity, defined through a different
// permittivity -- see function in_load(...). One can either use a mesh that is
// aligned to the load via curvilinear elements, or an unaligned mesh (parameter
// ALIGN_MESH). Convergence graphs saved both wrt. the dof number and cpu time.
//
// PDE: time-harmonic Maxwell's equations with discontinuous, piecewise-constant
//      electric permittivity
//      circular load in the middle of the large cavity, whose permittivity is
//      different from the rest of the domain
//
// Domain: square cavity with another small square cavity attached from outside
//         to its left edge
//
// Meshes: you can either use "oven_load_circle.mesh" containing curved elements
//         aligned with the circular load, or "oven_load_square.mesh" which is not
//         aligned.
//
// BC: perfect conductor on the boundary except for the right-most edge of the small
//     cavity, where a harmonic surface current is prescribed
//
// The following problem parameters can be changed easily:

const int P_INIT = 2;           // initial polynomial degree in mesh
const bool ALIGN_MESH = true;   // if ALIGN_MESH == true, curvilinear elements aligned with the
                                // circular load are used, otherwise one uses a non-aligned mesh.
const double THRESHOLD = 0.3;   // the adaptivity algorithm goes on until THRESHOLD*total_error is processed
                                // (see adapt_hcurl.cpp for explanation)
const int STRATEGY = 1;         // refinement strategy (0, 1, 2, 3 - see adapt_hcurl.cpp for explanation)
const bool H_ONLY = false;      // if H_ONLY == false then full hp-adaptivity takes place, otherwise
                                // h-adaptivity is used. Use this parameter to check that indeed adaptive
                                // hp-FEM converges much faster than adaptive h-FEM
const bool ISO_ONLY = false;    // when ISO_ONLY = true, only isotropic refinements are done,
                                // otherwise also anisotropic refinements are allowed
const int MESH_REGULARITY = -1; // specifies maximum allowed level of hanging nodes
                                // -1 ... arbitrary level hangning nodes
                                // 1, 2, 3,... means 1-irregular mesh, 2-irregular mesh, etc.
                                // total regularization (0) is not supported in adaptivity
const double ERR_STOP = 0.05;   // adaptivity process stops when error wrt. exact solution in H1 norm
                                // is less than this number
const int NDOF_STOP = 40000;    // adaptivity process stops when the number of degrees of freedom grows over
                                // this limit. This is mainly to prevent h-adaptivity to go on forever.

// remaining problem constants
const double e_0   = 8.8541878176 * 1e-12;
const double mu_0   = 1.256 * 1e-6;
const double e_r = 1.0;
const double mu_r = 1.0;
const double rho = 3820.0;
const double Cp = 7.531000;
const double freq = 1.0*2450000000.0;
const double omega = 2 * M_PI * freq;
const double c = 1 / sqrt(e_0 * mu_0);
const double kappa  = 2 * M_PI * freq * sqrt(e_0 * mu_0);
const double J = 0.0000033333;

int e_bc_types(int marker)
{
  if (marker == 2) return BC_ESSENTIAL; // perfect conductor
  else return BC_NATURAL; // impedance
}

bool in_load(double x, double y)
{
  double cx = -0.152994121;
  double cy =  0.030598824;
  double r = 0.043273273;
  if (sqr(cx - x) + sqr(cy - y) < sqr(r)) return true;
  else return false;
}

double gam(int marker, double x, double y)
{
  if (ALIGN_MESH && marker == 1) return 0.03;
  if (!ALIGN_MESH && in_load(x,y)) return 0.03;
  return 0.0;
}

double er(int marker, double x, double y)
{
  if (ALIGN_MESH && marker == 1) return 7.5;
  if (!ALIGN_MESH && in_load(x,y)) return 7.5;
  return 1.0;
}

void current1(double x, double y, scalar& F0, scalar& F1)
{
  F0 = 0.0;
  F1 = complex(0.0,  omega * J);
}

inline scalar int_fn_e_f(double (*fn)(int, double, double), RealFunction* fe, RealFunction* ff, RefMap* re, RefMap* rf)
{
  Quad2D* quad = re->get_quad_2d();
  int o = fe->get_fn_order() + ff->get_fn_order() + 2 + re->get_inv_ref_order();
  limit_order(o);
  fe->set_quad_order(o, FN_VAL);
  ff->set_quad_order(o, FN_VAL);

  int marker = re->get_active_element()->marker;
  double* x = re->get_phys_x(o);
  double* y = re->get_phys_y(o);

  double *e0 = fe->get_fn_values(0), *e1 = fe->get_fn_values(1);
  double *f0 = ff->get_fn_values(0), *f1 = ff->get_fn_values(1);

  scalar c;
  scalar result = 0.0;
  hcurl_integrate_jac_expression(( c = fn(marker, x[i], y[i]),
                                 c * (((*me)[0][0]*e0[i] + (*me)[0][1]*e1[i]) * ((*mf)[0][0]*f0[i] + (*mf)[0][1]*f1[i]) +
                                 ((*me)[1][0]*e0[i] + (*me)[1][1]*e1[i]) * ((*mf)[1][0]*f0[i] + (*mf)[1][1]*f1[i]))));
  return result;

}


// Maxwell
complex bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
   complex ikappa = complex(0.0, kappa);
   return   1.0/mu_r * int_curl_e_curl_f(fu, fv, ru, rv)
          - ikappa * sqrt(mu_0 / e_0) * int_fn_e_f(gam, fu, fv, ru, rv)
          - sqr(kappa) * int_fn_e_f(er, fu, fv, ru, rv);
}


complex linear_form_surf(RealFunction* fv, RefMap* refmap, EdgePos* ep)
{
  if (ep->marker == 2) return 0;
  return surf_int_F_f(current1, fv, refmap, ep);
}


int main(int argc, char* argv[])
{
  Mesh mesh;
  if (ALIGN_MESH) mesh.load("oven_load_circle.mesh");
  else mesh.load("oven_load_square.mesh");

  HcurlShapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  HcurlSpace space(&mesh, &shapeset);
  space.set_bc_types(e_bc_types);
  space.set_uniform_order(P_INIT);

  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form);
  wf.add_liform_surf(0, linear_form_surf);

  VectorView eview("Electric field",0,0,800, 590);
  OrderView ord("Order", 800, 0, 700, 590);

  UmfpackSolver umfpack;
  Solution sln, rsln;

  GnuplotGraph graph;
  graph.set_captions("Error Convergence for the Waveguide Problem", "Degrees of Freedom", "Error Estimate [%]");
  graph.add_row("error estimate", "-", "o");
  graph.set_log_y();

  GnuplotGraph graph_cpu;
  graph_cpu.set_captions("Error Convergence for the Waveguide Problem", "CPU Time", "Error Estimate [%]");
  graph_cpu.add_row("error estimate", "-", "o");
  graph_cpu.set_log_y();

  int it = 0;
  bool done = false;
  double cpu = 0.0;
  do
  {
    printf("\n---- it=%d ------------------------------------------------------------------\n\n", it++);
    begin_time();

    // enumerating basis functions
    space.assign_dofs();

    // coarse problem
    LinSystem sys(&wf, &umfpack);
    sys.set_spaces(1, &space);
    sys.set_pss(1, &pss);
    sys.assemble();
    sys.solve(1, &sln);

    cpu += end_time();

    // show real part of the solution
    AbsFilter abs(&sln);
    eview.set_min_max_range(0, 4e3);
    eview.show(&abs);
    ord.show(&space);

    // fine (reference) problem
    begin_time();
    RefSystem ref(&sys);
    ref.assemble();
    ref.solve(1, &rsln);

    // estimation of error and adaptation
    HcurlOrthoHP hp(1, &space);
    hp.set_kappa(sqr(kappa));
    double err_est = hp.calc_error(&sln, &rsln) * 100;
    info("Hcurl error estimate: %g%%", hcurl_error(&sln, &rsln) * 100);
    info("Adapt error estimate: %g%%", err_est);
    if (err_est < ERR_STOP || sys.get_num_dofs() >= NDOF_STOP) done = true;
    else hp.adapt(THRESHOLD, STRATEGY, H_ONLY, ISO_ONLY, MESH_REGULARITY);
    cpu += end_time();

    // plotting convergence wrt. numer of dofs
    graph.add_values(0, space.get_num_dofs(), err_est);
    graph.save("conv_dof.gp");

    // plotting convergence wrt. cpu time
    graph_cpu.add_values(0, cpu, err_est);
    graph_cpu.save("conv_cpu.gp");
  }
  while (done == false);
  verbose("\nTotal running time: %g sec", cpu);

  View::wait();
  return 0;
}

