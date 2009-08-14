#include "hermes2d.h"
#include "solver_umfpack.h"

// This example illustrates how to use nonhomogeneous (nonzero)
// Dirichlet boundary conditions.
//
// PDE: Poisson equation -Laplace u = CONST_F, where CONST_F is
// a constant right-hand side. It is not difficult to see that
// the function u(x,y) = (-CONST_F/4)*(x^2 + y^2) satisfies the
// above PDE. Since also the Dirichlet boundary conditions
// are chosen to match u(x,y), this function is the exact
// solution.
//
// Note that since the exact solution is a quadratic polynomial,
// Hermes will compute it exactly if all mesh elements are quadratic
// or higher (then the exact solution lies in the finite element space).
// If some elements in the mesh are linear, Hermes will only find
// an approximation, Below you can play with the parameters CONST_F,
// P_INIT, and UNIFORM_REF_LEVEL.

double CONST_F = -4.0;       // constant right-hand side
int P_INIT = 2;              // initial polynomial degree in all elements
int UNIFORM_REF_LEVEL = 3;   // number of initial uniform mesh refinements

// boundary condition types
// Note: essential means Dirichlet (prescribed is value of solution at the boundary).
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// function values for Dirichlet boundary markers
scalar bc_values(int marker, double x, double y)
{
  return (-CONST_F/4.0)*(x*x + y*y);
}

// bilinear form
scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}

// linear form
scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return CONST_F*int_v(fv, rv);
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

  // wait for keyboard or mouse input
  printf("Waiting for keyboard or mouse input.\n");
  View::wait();
  return 0;
}
