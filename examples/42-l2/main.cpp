#include "hermes2d.h"
#include "solver_umfpack.h"

//  This example tests L2 space and L2 shapeset.
//  It projects continuous function x^3 + y^3 onto L2 space in L2 norm.
//  When zero order used, solution is a piecewice constant function.

scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_u_v(fu, fv, ru, rv);
}

double F(double x, double y)
{
  return x*x*x + y*y*y;
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return int_F_v(F, fv, rv);
}

int bc_types(int marker)
{
   return BC_NONE;
}

int P_INIT = 0;

int main(int argc, char* argv[])
{
  if (argc < 2) error("Missing mesh parameter.");

  Mesh mesh;
  mesh.load(argv[1]);
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  L2Shapeset shapeset;
  L2Space space(&mesh, &shapeset);
  PrecalcShapeset pss(&shapeset);

  space.set_bc_types(bc_types);
  space.set_uniform_order(P_INIT);
  space.assign_dofs();

/*  BaseView bview;
  bview.show(&space);
  bview.wait_for_close();*/

  Solution sln;
  UmfpackSolver umfpack;

  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form);
  wf.add_liform(0, linear_form);

  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);

  sys.assemble();
  sys.solve(1, &sln);

  ScalarView view1("Solution 1");
  view1.show(&sln);
  view1.wait_for_keypress();

  View::wait();
  return 0;
}

