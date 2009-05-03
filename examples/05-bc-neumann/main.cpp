#include "hermes2d.h"
#include "solver_umfpack.h"

// This example is a continuation of examples 03 and 04. It shows 
// you how to use Neumann boundary conditions.
//
// PDE: Poisson equation -Laplace u = f, where f = -1 
//
// The exact solution to this problem has a singular (infinite) 
// gradient at the re-entrant corner. You can play with the 
// initial polynomial degree P_INIT, with the number of refinements 
// towards the re-entrant corner, etc. 

int P_INIT = 2;             // initial polynomial degree in all elements
int CORNER_REF_LEVEL = 10;  // number of mesh refinements towards the re-entrant corner

int bc_types(int marker)
{
  // Note: essential means Dirichlet (prescribed is value of solution at the boundary). 
  // Natural means Neumann, Newton, or any other type of condition where the solution 
  // value is not prescribed.
  return (marker == 4) ? BC_ESSENTIAL : BC_NATURAL; 
}

scalar bc_values(int marker, double x, double y)
{
  if (marker == 4) return 0.0; // Dirichlet BC value
  else return 0.0;             // Neumann BC value
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
  mesh.refine_towards_vertex(3, CORNER_REF_LEVEL);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
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

  // visualize the approximation
  ScalarView view("Solution", 0, 0, 600, 600);
  view.show(&sln);

  // TODO: show gradient magnitude



 
  printf("Waiting for keyboard or mouse input.\n");
  View::wait();
  return 0;
}
