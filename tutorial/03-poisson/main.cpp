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

int P_INIT = 5;            // Uniform polynomial degree of mesh elements.

// Problem parameters.
double CONST_F = 2.0;  

// Boundary condition types.
// Note: "essential" boundary condition means that 
// the solution value is prescribed.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements.
  mesh.refine_element(0);

  // Initialize the shapeset.
  H1Shapeset shapeset;

  // Create an H1 space.
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_essential_bc_values(essential_bc_values);
  space.set_uniform_order(P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_biform(callback(bilinear_form));
  wf.add_liform(callback(linear_form));

  // Matrix solver.
  UmfpackSolver solver;

  // Initialize the linear system.
  LinSystem ls(&wf, &solver, &space);

  // Assemble and solve the matrix problem.
  Solution sln;
  ls.assemble();
  ls.solve(&sln);

  // Visualize the solution.
  ScalarView view("Solution");
  view.show(&sln);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

