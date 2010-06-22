#include "hermes2d.h"

// This example shows how to define Neumann boundary conditions. In addition,
// you will see how a Filter is used to visualize gradient of the solution
//
// PDE: Poisson equation -Laplace u = f, where f = CONST_F
//
// BC: u = 0 on Gamma_4 (edges meeting at the re-entrant corner)
//     du/dn = CONST_GAMMA_1 on Gamma_1 (bottom edge)
//     du/dn = CONST_GAMMA_2 on Gamma_2 (top edge, circular arc, and right-most edge)
//     du/dn = CONST_GAMMA_3 on Gamma_3 (left-most edge)
//
// You can play with the parameters below. For most choices of the four constants,
// the solution has a singular (infinite) gradient at the re-entrant corner.
// Therefore we visualize not only the solution but also its gradient.

int P_INIT = 4;               // Initial polynomial degree in all elements.
int CORNER_REF_LEVEL = 12;    // Number of mesh refinements towards the re-entrant corner.

// Problem parameters.
double CONST_F = -1.0;                        // Right-hand side.
double CONST_GAMMA[3] = {-0.5, 1.0, -0.5};    // Outer normal derivative on Gamma_1,2,3.

// Boundary condition types.
// Note: "natural" boundary condition means that 
// the solution value is not prescribed.
BCType bc_types(int marker)
{
  return (marker == 4) ? BC_ESSENTIAL : BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0.0;
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
  mesh.refine_towards_vertex(3, CORNER_REF_LEVEL);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));
  wf.add_vector_form_surf(callback(linear_form_surf));

  // Initialize the linear system.
  LinSystem ls(&wf, &space);

  // Assemble and solve the matrix problem.
  Solution sln;
  ls.assemble();
  ls.solve(&sln);

  // Visualize the approximation.
  ScalarView view("Solution", 0, 0, 600, 600);
  view.show(&sln);

  // Compute and show gradient magnitude.
  // (Note that the gradient at the re-entrant
  // corner needs to be truncated for visualization purposes.)
  ScalarView gradview("Gradient", 650, 0, 600, 600);
  MagFilter grad(&sln, &sln, H2D_FN_DX, H2D_FN_DY);
  gradview.show(&grad);

  // Wait for the views to be closed.
  View::wait();
  return 0;
}
