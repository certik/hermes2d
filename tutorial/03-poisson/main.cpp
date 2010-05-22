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

// boundary condition types (essential = Dirichlet)
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// function values for essential(Dirichlet) boundary conditions
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

// return the value \int \nabla u . \nabla v dx
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

// return the value \int v dx
template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return CONST_F * int_v<Real, Scalar>(n, wt, v);
}

int main(int argc, char* argv[])
{
  // Load the mesh file.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform sample initial mesh refinement.
  mesh.refine_element(0);

  // Initialize the shapeset and the cache.
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // Create an H1 space.
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_essential_bc_values(essential_bc_values);
  space.set_uniform_order(P_INIT);

  // Enumerate degrees of freedom.
  int ndof = assign_dofs(&space);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_biform(callback(bilinear_form));
  wf.add_liform(callback(linear_form));

  // Initialize the linear system and solver.
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_space(&space);
  sys.set_pss(&pss);

  // Assemble the stiffness matrix and solve the system.
  Solution sln;
  sys.assemble();
  sys.solve(&sln);

  // Visualize the solution.
  ScalarView view("Solution");
  view.show(&sln);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

