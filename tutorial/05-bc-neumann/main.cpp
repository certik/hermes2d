#include "hermes2d.h"
#include "solver_umfpack.h"

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

double CONST_F = -1.0;        // right-hand side
double CONST_GAMMA[3] = {-0.5, 1.0, -0.5}; // outer normal derivative on Gamma_1,2,3
int P_INIT = 4;               // initial polynomial degree in all elements
int CORNER_REF_LEVEL = 12;    // number of mesh refinements towards the re-entrant corner

// boundary condition types
// Note: natural means Neumann, Newton, or any other type of condition
// where the solution value is not prescribed.
int bc_types(int marker)
{
  return (marker == 4) ? BC_ESSENTIAL : BC_NATURAL;
}

// function values for Dirichlet boundary markers
scalar bc_values(int marker, double x, double y)
{
  return 0.0;
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return CONST_F*int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return CONST_GAMMA[e->marker - 1] * int_v<Real, Scalar>(n, wt, v);
}

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);
  mesh.refine_towards_vertex(3, CORNER_REF_LEVEL);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
  
  // enumerate degrees of freedom
  int ndof = assign_dofs(&space);

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form));
  wf.add_liform(0, callback(linear_form));
  wf.add_liform_surf(0, callback(linear_form_surf));

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

  // compute and show gradient magnitude
  // (note that the infinite gradient at the re-entrant
  // corner needs to be truncated for visualization purposes)
  ScalarView gradview("Gradient", 650, 0, 600, 600);
  MagFilter grad(&sln, &sln, H2D_FN_DX, H2D_FN_DY);
  gradview.show(&grad);

  // wait for a view to be closed
  View::wait();
  return 0;
}
