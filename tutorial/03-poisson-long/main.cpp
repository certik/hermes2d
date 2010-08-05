#define H2D_REPORT_INFO
#include "hermes2d.h"

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

int P_INIT = 3;                                   // Uniform polynomial degree of mesh elements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
double CONST_F = 2.0;  

// Boundary condition types.
// Note: "essential" means that solution value is prescribed.
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

  // Perform initial mesh refinements (optional).
  //mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));

  // Initialize the linear problem.
  LinearProblem lp(&wf, &space);

  // Select matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;
  init_matrix_solver(matrix_solver, ndof, mat, rhs, solver);

  // Assemble stiffness matrix and rhs.
  lp.assemble(mat, rhs);

  // Solve the matrix problem.
  if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

  // Convert coefficient vector into a Solution.
  Solution sln;
  sln.set_fe_solution(&space, rhs);

  // Visualize the solution.
  WinGeom* sln_win_geom = new WinGeom(0, 0, 440, 350);
  ScalarView view("Solution", sln_win_geom);
  view.show(&sln);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

