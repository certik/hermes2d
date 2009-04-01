#include "hermes2d.h"
#include "solver_umfpack.h"


scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return int_v(fv, rv);
}

void der(int n, scalar* a, scalar* out)
{
  for (int i = 0; i < n; i++)
  {
    out[i] = a[i];
  }
}

int main(int argc, char* argv[])
{
  Mesh mesh;
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  H1Space space(&mesh, &shapeset);
  Solution sln;
  UmfpackSolver umfpack;

  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form);
  wf.add_liform(0, linear_form);

  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);

  mesh.load("square1.mesh");
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  space.set_uniform_order(3);
  space.assign_dofs();
 
  sys.assemble();
  sys.solve(1, &sln);

  ScalarView view1("Solution 1");
  view1.show(&sln);
  view1.wait_for_keypress();

  ////////////////////////

  mesh.load("square2.mesh");
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  
  space.set_uniform_order(10);
  space.assign_dofs();

  sys.assemble();
  sys.solve(1, &sln);

  ScalarView view2("Solution 2");
  view2.show(&sln);

  ////////////////////////
  
  mesh.load("diamond.mesh");
  mesh.refine_all_elements(1);
  mesh.refine_all_elements(1);
  
  space.set_uniform_order(10);
  space.assign_dofs();

  sys.assemble();
  sys.solve(1, &sln);

  ScalarView view3("Solution 3");
  view3.show(&sln);

  View::wait();
  return 0;
}
