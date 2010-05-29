#include "hermes2d.h"
#include "solver_umfpack.h"

// This example explains how to create two spaces over a mesh and use them
// to solve a simple problem of linear elasticity. At the end, VonMises
// filter is used to visualize the stress.
//
// PDE: Lame equations of linear elasticity.
//
// BC: du_1/dn = f_0 on Gamma_3 and du_1/dn = 0 on Gamma_2, Gamma_4, Gamma_5,
//     du_2/dn = f_1 on Gamma_3 and du_2/dn = 0 on Gamma_2, Gamma_4, Gamma_5,
//     u_1 = 0 and u_2 = 0 on Gamma_1.
//
// The following parameters can be changed:

const int P_INIT = 8;                                      // Initial polynomial degree of all elements.

// Problem parameters.
const double E  = 200e9;                                   // Young modulus (steel).
const double nu = 0.3;                                     // Poisson ratio.
const double f_0  = 0;                                     // External force in x-direction.
const double f_1  = 1e4;                                   // External force in y-direction.
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));  // First Lame constant.
const double mu = E / (2*(1 + nu));                        // Second Lame constant.

// Boundary marker (external force).
const int GAMMA_3_BDY = 3;

// Boundary condition types.
BCType bc_types(int marker)
  { return (marker == 1) ? BC_ESSENTIAL : BC_NATURAL; }

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
  { return 0; }

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("sample.mesh", &mesh);

  // Initialize the shapeset and the cache.
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // Create the x displacement space.
  H1Space xdisp(&mesh, &shapeset);
  xdisp.set_bc_types(bc_types);
  xdisp.set_essential_bc_values(essential_bc_values);
  xdisp.set_uniform_order(P_INIT);

  // Create the y displacement space.
  H1Space ydisp(&mesh, &shapeset);
  ydisp.set_bc_types(bc_types);
  ydisp.set_essential_bc_values(essential_bc_values);
  ydisp.set_uniform_order(P_INIT);

  // Enumerate degrees of freedom.
  int ndof = assign_dofs(2, &xdisp, &ydisp);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_biform(0, 0, callback(bilinear_form_0_0), H2D_SYM);  // Note that only one symmetric part is
  wf.add_biform(0, 1, callback(bilinear_form_0_1), H2D_SYM);  // added in the case of symmetric bilinear
  wf.add_biform(1, 1, callback(bilinear_form_1_1), H2D_SYM);  // forms.
  wf.add_liform_surf(0, callback(linear_form_surf_0), GAMMA_3_BDY);
  wf.add_liform_surf(1, callback(linear_form_surf_1), GAMMA_3_BDY);

  // Matrix solver.
  UmfpackSolver umfpack;

  // Initialize the linear system.
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(2, &xdisp, &ydisp);
  sys.set_pss(&pss);

  // Assemble and solve the matrix problem.
  Solution xsln, ysln;
  sys.assemble();
  sys.solve(2, &xsln, &ysln);

  // Visualize the solution.
  ScalarView view("Von Mises stress [Pa]", 50, 50, 1200, 600);
  VonMisesFilter stress(&xsln, &ysln, lambda, mu);
  view.show(&stress, H2D_EPS_HIGH, H2D_FN_VAL_0, &xsln, &ysln, 1.5e5);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

