#include "hermes2d.h"
#include "solver_umfpack.h"  // defines the class UmfpackSolver


scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  // return the value \int \nabla u . \nabla v dx
  return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  // return the value \int v dx
  return 2*int_v(fv, rv);
}


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");
  mesh.refine_element(0);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_uniform_order(5);
  space.assign_dofs();

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form);
  wf.add_liform(0, linear_form);

  // initialize the linear system and solver
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);
  
  // assemble the stiffness matrix and solve the system
  Solution sln;
  sys.assemble();
  sys.solve(1, &sln);

  // visualize the solution
  ScalarView view("Solution");
  view.show(&sln);
  
  View::wait();
  return 0;
}
