#include "hermes2d.h"
#include "solver_umfpack.h"


int bc_types(int marker)
{
  return BC_NATURAL;
}

scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return -4*int_v(fv, rv);
}

scalar linear_form_surf(RealFunction* fv, RefMap* rv, EdgePos* ep)
{
  return 2*surf_int_v(fv, rv, ep);
}


int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("square.mesh");

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_uniform_order(5);
  space.fix_vertex(0, 1.0);
  space.assign_dofs();

  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form);
  wf.add_liform(0, linear_form);
  wf.add_liform_surf(0, linear_form_surf, ANY);
  
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);

  Solution sln;
  sys.assemble();
  sys.solve(1, &sln);
  
  ScalarView view("Solution");
  view.show(&sln);
  
  View::wait();
  return 0;
}

