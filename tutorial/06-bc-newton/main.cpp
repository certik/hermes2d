#include "hermes2d.h"
#include "solver_umfpack.h"

// This example explains how to use Newton boundary conditions. Again,
// a Filter is used to visualize the solution gradient.
//
// PDE: Laplace equation -Laplace u = 0 (this equation describes, among
// many other things, also stationary heat transfer in a homogeneous linear
// material).
//
// BC: u = T1 ... fixed temperature on Gamma_3 (Dirichlet)
//     du/dn = 0 ... insulated wall on Gamma_2 and Gamma_4 (Neumann)
//     du/dn = H*(u - T0) ... heat flux on Gamma_1 (Newton)
//
// Note that the last BC can be written in the form  du/dn - H*u = -h*T0.
//
// You can play with the parameters below:

double T1 = 30.0;            // prescribed temperature on Gamma_3
double T0 = 20.0;            // outer temperature on Gamma_1
double H  = 0.05;            // heat flux on Gamma_1
int P_INIT = 6;              // uniform polynomial degree in the mesh
int UNIFORM_REF_LEVEL = 2;   // number of initial uniform mesh refinements
int CORNER_REF_LEVEL = 12;   // number of mesh refinements towards the re-entrant corner

// boundary condition types
int bc_types(int marker)
  { return (marker == 3) ? BC_ESSENTIAL : BC_NATURAL; }

// function values for Dirichlet boundary markers
scalar bc_values(int marker, double x, double y)
  { return T1; }

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return H * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return T0 * H * int_v<Real, Scalar>(n, wt, v);
}


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);
  for(int i=0; i<UNIFORM_REF_LEVEL; i++) mesh.refine_all_elements();
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
  wf.add_biform_surf(0, 0, callback(bilinear_form_surf), 1);
  wf.add_liform_surf(0, callback(linear_form_surf), 1);

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
  ScalarView view("Solution", 0, 0, 600, 600);
  view.show(&sln);

  // compute and show gradient magnitude
  // (note that the infinite gradient at the re-entrant
  // corner needs to be truncated for visualization purposes)
  ScalarView gradview("Gradient", 650, 0, 600, 600);
  MagFilter grad(&sln, &sln, H2D_FN_DX, H2D_FN_DY);
  gradview.show(&grad);

  // wait for all views to be closed
  View::wait();
  return 0;
}
