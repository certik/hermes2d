#include "hermes2d.h"
#include "solver_umfpack.h"

// This example explains how to create two spaces over a mesh and use them
// to solve a simple problem of linear elasticity. At the end, VonMises
// filter is used to visualize the stress.
//
// PDE: Lame equations of linear elasticity
//
// BC: du_1/dn = f_0 on Gamma_3 and du_1/dn = 0 on Gamma_2, Gamma_4, Gamma_5
//     du_2/dn = f_1 on Gamma_3 and du_2/dn = 0 on Gamma_2, Gamma_4, Gamma_5
//     u_1 = 0 and u_2 = 0 on Gamma_1
//
// The following parameters can be changed:

const int P_INIT = 8;                                      // initial polynomial degree in all elements

// problem constants
const double E  = 200e9;                                   // Young modulus (steel)
const double nu = 0.3;                                     // Poisson ratio
const double f_0  = 0;                                     // external force in x-direction
const double f_1  = 1e4;                                   // external force in y-direction
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));  // first Lame constant
const double mu = E / (2*(1 + nu));                        // second Lame constant

// boundary marker for the external force
const int GAMMA_3_BDY = 3;

// boundary condition types
int bc_types(int marker)
  { return (marker == 1) ? BC_ESSENTIAL : BC_NATURAL; }

// function values for Dirichlet boundary conditions
double bc_values(int marker, double x, double y)
  { return 0; }

// bilinear forms
template<typename Real, typename Scalar>
Scalar bilinear_form_0_0(int n, double *wt, Func<Real> *u, Func<Real> *v,
                         Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (lambda + 2*mu) * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
                      mu * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_0_1(int n, double *wt, Func<Real> *u, Func<Real> *v,
                         Geom<Real> *e, ExtData<Scalar> *ext)
{
  return lambda * int_dudy_dvdx<Real, Scalar>(n, wt, u, v) +
             mu * int_dudx_dvdy<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_1(int n, double *wt, Func<Real> *u, Func<Real> *v,
                         Geom<Real> *e, ExtData<Scalar> *ext)
{
  return              mu * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
         (lambda + 2*mu) * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
}

// linear forms
template<typename Real, typename Scalar>
Scalar linear_form_surf_0(int n, double *wt, Func<Real> *v, Geom<Real> *e,
                          ExtData<Scalar> *ext)
{
  return f_0 * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_1(int n, double *wt, Func<Real> *v, Geom<Real> *e,
                          ExtData<Scalar> *ext)
{
  return f_1 * int_v<Real, Scalar>(n, wt, v);
}

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("sample.mesh", &mesh);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create the x displacement space
  H1Space xdisp(&mesh, &shapeset);
  xdisp.set_bc_types(bc_types);
  xdisp.set_bc_values(bc_values);
  xdisp.set_uniform_order(P_INIT);

   // create the y displacement space
  H1Space ydisp(&mesh, &shapeset);
  ydisp.set_bc_types(bc_types);
  ydisp.set_bc_values(bc_values);
  ydisp.set_uniform_order(P_INIT);

  // enumerate degrees of freedom
  int ndof = assign_dofs(2, &xdisp, &ydisp);

  // initialize the weak formulation
  WeakForm wf(2);
  wf.add_biform(0, 0, callback(bilinear_form_0_0), H2D_SYM);  // Note that only one symmetric part is
  wf.add_biform(0, 1, callback(bilinear_form_0_1), H2D_SYM);  // added in the case of symmetric bilinear
  wf.add_biform(1, 1, callback(bilinear_form_1_1), H2D_SYM);  // forms.
  wf.add_liform_surf(0, callback(linear_form_surf_0), GAMMA_3_BDY);
  wf.add_liform_surf(1, callback(linear_form_surf_1), GAMMA_3_BDY);

  // initialize the linear system and solver
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(2, &xdisp, &ydisp);
  sys.set_pss(1, &pss);

  // assemble the stiffness matrix and solve the system
  Solution xsln, ysln;
  sys.assemble();
  sys.solve(2, &xsln, &ysln);

  // visualize the solution
  ScalarView view("Von Mises stress [Pa]", 50, 50, 1200, 600);
  VonMisesFilter stress(&xsln, &ysln, lambda, mu);
  view.show(&stress, H2D_EPS_HIGH, H2D_FN_VAL_0, &xsln, &ysln, 1.5e5);

  // wait for a view to be closed
  View::wait();
  return 0;
}

