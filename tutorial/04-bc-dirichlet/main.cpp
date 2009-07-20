#include "hermes2d.h"
#include "solver_umfpack.h"

// This example is a continuation of example 03. It shows you how
// to use nonhomogeneous (nonzero) Dirichlet boundary conditions.
//
// PDE: Poisson equation -Laplace u = f, where f = -4 corresponds
// to the function x^2 + y^2. Since also the Dirichlet boundary
// conditions correspond to the same function, u(x) = x^2 + y^2
// is the exact solution of this problem.
//
// Note that since the exact solution is a quadratic polynomial,
// Hermes will compute it exactly if all mesh elements have polynomial
// degree at least 2 (because then the exact solution lies in the
// finite element space). If you choose at least one element to be
// linear, Hermes will only find an approximation, You can try this
// easily, just redefine below P_INIT to 1. You can also play with
// the number of initial uniform refinements UNIFORM_REF_LEVEL.

int P_INIT = 2;              // initial polynomial degree in all elements
int UNIFORM_REF_LEVEL = 3;   // number of initial uniform mesh refinements

int bc_types(int marker)
{
  // all markers denote the essential (Dirichlet) boundary condition
  return BC_ESSENTIAL;
}

scalar bc_values(int marker, double x, double y)
{
  // this is the Dirichlet BC value for all markers
  return x*x + y*y;
}

scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return -4*int_v(fv, rv);
}


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");
  for(int i=0; i<UNIFORM_REF_LEVEL; i++) mesh.refine_all_elements();

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

  printf("Waiting for keyboard or mouse input.\n");
  View::wait();
  return 0;
}
