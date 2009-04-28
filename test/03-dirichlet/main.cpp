#include "hermes2d.h"
#include "solver_umfpack.h"


int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

double bc_values(int marker, double x, double y)
{
  if (marker == 1) return 1.0 - sqr(x);
  return 0.0;
}


scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}

double xx = 1.0;

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return xx * int_v(fv, rv);
}


int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("square2.mesh");
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(5);
  space.assign_dofs();

  BaseView bview;
  bview.show(&space);

  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form);
  wf.add_liform(0, linear_form);

  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);

  Solution sln;
  sys.assemble();
  sys.solve(1, &sln);

  ScalarView view("Solution", 200, 150, 1000, 900);
  view.show(&sln);
  
  /*view.wait_for_keypress();
  sln.multiply(2);
  view.show(&sln);*/
  
  /*view.wait_for_keypress();
  xx = 3.0;
  sys.assemble_rhs_only();
  sys.solve(1, &sln);
  view.show(&sln);*/

  View::wait();
  return 0;
}
