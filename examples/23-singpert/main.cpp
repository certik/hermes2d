#include "hermes2d.h"
#include "solver_umfpack.h"

//
//  PDE: -Laplace u + K*K*u = const
//
//  Domain: square, see the file singpert.mesh
//
//  BC:  Homogeneous Dirichlet
//
//  With large K, this is a singularly perturbed problem that exhibits an extremely
//  thin and steep boundary layer. Singularly perturbed problems are considered to
//  be very difficult, but you'll see that Hermes performs not bad at all. We tried
//  various values of K between 1 and 10^10, always starting from an extremely coarse
//  mesh. The only problem that we noticed was that for large K, all refinements go
//  towards the edges for a long time, and one has to wait a lot before the mesh gets
//  refined in the corners. To keep the solution reasonably scaled, set F_CONST to be
//  roughly K*K. During each computation, an approximate convergence curve is saved
//  to the file "convergence.gp". As in other adaptivity examples, you can compare
//  hp-adaptivity (H_ONLY=0) with h-adaptivity (H_ONLY=1). Moreover, you can turn off
//  and on anisotropic element refinements via the ISO_ONLY parameter, and see how
//  important they are for the efficient resolution of boundary layers. For large K
//  and with ISO_ONLY = 0, you'll see many levels of hanging nodes. The following
//  problem parameters can be changed easily:
//

double K = 1e3;           // equation parameter
double CONST_F = 1e6;     // constant right-hand side (set to be roughly K*K for scaling purposes)
int P_INIT = 1;           // initial polynomial degree in mesh
int INIT_REF_NUM = 1;     // number of initial mesh refinements (the original mesh is just one element)
int INIT_REF_NUM_BDY = 0; // number of initial mesh refinements towards the boundary
double THRESHOLD = 0.3;   // the adaptivity algorithm goes on until THRESHOLD*total_error is processed
                          // (see adapt_h1.cpp for explanation)
int STRATEGY = 0;         // refinement strategy (0, 1, 2, 3 - see adapt_h1.cpp for explanation)
int H_ONLY = 0;           // if H_ONLY == 0 then full hp-adaptivity takes place, otherwise
                          // h-adaptivity is used. Use this parameter to check that indeed adaptive
                          // hp-FEM converges much faster than adaptive h-FEM
int ISO_ONLY = 0;         // with ISO_ONLY = 0, quadrilateral elements can be refined isotropically or
                          // anisotropically, otherwise they only can be refined isotropically. Check
                          // that the adaptivity performs much better with ISO_ONLY = 0!
double ERR_STOP = 0.01;   // adaptivity process stops when error wrt. exact solution in H1 norm
                          // is less than this number
int NDOF_STOP = 100000;   // adaptivity process stops when the number of degrees of freedom grows over
                          // this limit. This is mainly to prevent h-adaptivity to go on forever.

scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
   // return the value \int \nabla u . \nabla v dx
  return int_grad_u_grad_v(fu, fv, ru, rv) + K*K*int_u_v(fu, fv, ru, rv);

}

double rhs(double x, double y)
{
  return CONST_F;
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return int_F_v(rhs, fv, rv);
}

int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("singpert.mesh");

  // initial mesh refinement (here you can apply arbitrary
  // other initial refinements, see example 01)
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, INIT_REF_NUM_BDY);

  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
  //space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form, SYM);
  wf.add_liform(0, linear_form);

  ScalarView sview("Coarse solution", 0, 100, 798, 700);
  OrderView  oview("Polynomial orders", 800, 100, 798, 700);

  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Error Convergence for the Singularly Perturbed Problem", "Degrees of Freedom", "Error [%]");
  graph.add_row("error estimate", "k", "--");

  Solution sln, rsln;
  UmfpackSolver umfpack;

  int it = 1;
  begin_time();
  while (1)
  {
    info("\n---- it=%d ------------------------------------------------------------------\n", it++);

    space.assign_dofs();

    LinSystem ls(&wf, &umfpack);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln);

    sview.show(&sln);
    oview.show(&space);

    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &rsln);

    H1OrthoHP hp(1, &space);
    double err_est = hp.calc_error(&sln, &rsln) * 100;
    info("Estimate of error: %g%%", err_est);

    graph.add_values(0, space.get_num_dofs(), err_est);
    graph.save("convergence.gp");

    if (err_est < ERR_STOP || ls.get_num_dofs() >= NDOF_STOP) break;
    hp.adapt(THRESHOLD, STRATEGY, H_ONLY, ISO_ONLY);
  }

  verbose("\nTotal run time: %g sec", end_time());
  View::wait();
  return 0;
}
