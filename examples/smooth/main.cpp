#include "hermes2d.h"
#include "solver_umfpack.h"

//
//  This example shows that the worst thing you can ever do is to approximate
//  smooth parts of solutions with uniform low-order meshes. The exact solution
//  to this Poisson problem is u(x,y) = sin(x)*sin(y), defined in the square
//  (0, pi)x(0, pi). Below, set H_REFIN = false for pure p-adaptivity (without
//  any spatial refinements), and H_REFIN = true for h-adaptivity. Compare the
//  convergence curves.
//
//  PDE: -Laplace u = f
//
//  Known exact solution, see functions fn() and fndd()
//
//  Domain: square domain (0, pi)x(0, pi), mesh files square_quad.mesh, square_tri.mesh
//
//  BC:  homogeneous Dirichlet
//
//  TODO: Implement preconditioned CG. Then it should be possible to go much higher in the
//  number of DOF for the low-order approximation.


const int P_INIT = 1;             // initial polynomial degree in mesh
const bool H_REFIN = false;        // if H_REFIN == true then Hermes will do uniform h-refinements,
                                  // otherwise p-refinements
const int INIT_REF_NUM = 1;       // number of initial uniform mesh refinements: use 0 for one
                                  // element only, 1 for 4 elements, 2 for 15 elements, etc.

// exact solution (for error calculation only)
static double fn(double x, double y)
{
  return sin(x)*sin(y);
}

// exact solution derivatives (for error calculation only)
static double fndd(double x, double y, double& dx, double& dy)
{
  dx = cos(x)*sin(y);
  dy = sin(x)*cos(y);
  return fn(x, y);
}

scalar bc_values(int marker, double x, double y)
{
  return 0;
}

scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}

double rhs(double x, double y)
{
  return 2*sin(x)*sin(y);
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return int_F_v(rhs, fv, rv);
}

int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("square_quad.mesh");
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  if(P_INIT < 2 && INIT_REF_NUM < 1) {
    info("No degrees of freedom in the discrete problem.\nEither subdivide the mesh or increase p.\nExiting.\n");
    exit(0);
  }

  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
  space.set_bc_values(bc_values);
  int actual_poly_degree = P_INIT;
  space.set_uniform_order(actual_poly_degree);

  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form, SYM);
  wf.add_liform(0, linear_form);

  ScalarView sview("Coarse solution", 0, 100, 798, 700);
  OrderView  oview("Polynomial orders", 800, 100, 798, 700);

  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Error Convergence for the Smooth Problem", "Degrees of Freedom", "Error [%]");
  graph.add_row("exact error", "k", "-", "o");

  Solution sln, rsln;
  UmfpackSolver umfpack;

  int it = 1;
  bool done = false;
  do
  {
    info("\n---- it=%d ------------------------------------------------------------------\n", it++);

    // enumerating basis functions
    space.assign_dofs();

    LinSystem ls(&wf, &umfpack);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln);

    // calculate error wrt. exact solution
    ExactSolution exact(&mesh, fndd);
    double error = h1_error(&sln, &exact) * 100;
    info("Exact solution error: %g%%", error);

    // plotting convergence wrt. number of dofs
    graph.add_values(0, space.get_num_dofs(), error);
    graph.save("conv_dof.gp");

    // view the solution
    sview.show(&sln);
    oview.show(&space);
    info("Click into the solution window and press spacebar to continue.");
    sview.wait_for_keypress();

    // refine the mesh uniformly either in 'h' or 'p'
    if (H_REFIN == true) {
      mesh.refine_all_elements();
      space.set_uniform_order(actual_poly_degree);
    }
    else {
      actual_poly_degree++;
      if(actual_poly_degree > 10) done = true;
      else space.set_uniform_order(actual_poly_degree);
    }
  }
  while (done == false);

  // wait for keyboard or mouse input
  printf("Waiting for keyboard or mouse input.\n");
  View::wait();
  return 0;
}
