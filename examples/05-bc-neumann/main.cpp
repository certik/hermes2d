#include "hermes2d.h"
#include "solver_umfpack.h"


int bc_types(int marker)
{
  return (marker == 4) ? BC_ESSENTIAL : BC_NATURAL;
}

scalar bc_values(int marker, double x, double y)
{
  return 0.0; // Dirichlet BC value
}


scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return -int_v(fv, rv);
}

scalar linear_form_surf_Gamma_2(RealFunction* fv, RefMap* rv, EdgePos* ep)
{
  return 1.0 * surf_int_v(fv, rv, ep);
}


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(5);
  space.assign_dofs();

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form);
  wf.add_liform(0, linear_form);
  wf.add_liform_surf(0, linear_form_surf_Gamma_2, 2);
  
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
