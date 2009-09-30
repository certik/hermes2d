#include "hermes2d.h"
#include "solver_umfpack.h"  // defines the class UmfpackSolver

// This example shows how to solve a first simple PDE:
//   - load the mesh,
//   - perform initial refinements
//   - create a H1 space over the mesh
//   - define weak formulation
//   - initialize matrix solver
//   - assemble and solve the matrix system
//   - visualize the solution
//
// PDE: Poisson equation -Laplace u = CONST_F with homogeneous (zero)
//      Dirichlet boundary conditions.
//
// You can change the constant right-hand side CONST_F, the
// initial polynomial degree P_INIT, and play with various initial
// mesh refinements at the beginning of the main() function.

double CONST_F = 2.0;   // Constant right-hand side.
int P_INIT = 5;         // Uniform polynomial degree of mesh elements.

// bilinear form
scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv); // returns \int \nabla u . \nabla v dx
}

// linear form
scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return CONST_F*int_v(fv, rv); // returns \int v dx
}

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");

  // sample "manual" mesh refinement
  mesh.refine_element(0);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_uniform_order(P_INIT);
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

  // wait for keyboard or mouse input
  printf("Waiting for keyboard or mouse input.\n");
  View::wait();
  return 0;
}
