#include "hermes2d.h"
#include "solver_umfpack.h"

// This example shows how a multigroup neutron diffusion equation with 4 groups in the reactor core
// can be solved using Hermes2d. Eigenproblem is solved using power interations,
// in each iteration a PDE is solved using FEM.
//
// The reactor neutronics is given by the following eigenproblem:
//
//  - \nabla \cdot D_g \nabla \phi_g + \Sigma_{Rg}\phi_g - \sum_{g' \neq g} \Sigma_s^{g'\to g} \phi_{g'} =
//  = \frac{\chi_g}{k_{eff}} \sum_{g'} \nu_{g'} \Sigma_{fg'}\phi_{g'}
//
// where 1/k_{eff} is eigenvalue and \phi_g, g = 1,...,4 are eigenvectors (neutron fluxes)
//
// BC:
//
// homogeneous neumann on symmetry axis
// d \phi_g / d n = - 0.5 \phi_g   elsewhere
//
// The eigenproblem is numerically solved using common technique known as the power method (power iterations):
//
//  1) Make an initial estimate of \phi_g and k_{eff}
//  2) For n = 1, 2,...
//         solve for \phi_g using previous k_prev
//         solve for new k_{eff}
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{prev}
//               k_new =  k_prev -------------------------------------------------------------------------
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{new}
//  3) Stop iterations when
//
//     |   k_new - k_prev  |
//     | ----------------- |  < epsilon
//     |       k_new       |
//
//


const int init_order = 1;
const int init_refinements = 4;

// Boundary and area markers
const int marker_reflector = 1;
const int marker_core = 2;

const int bc_vacuum = 1;
const int bc_sym = 2;

// Boundary condition
int bc_types(int marker)
{
  return BC_NATURAL;
}

// reflector properties (0) core properties (1),
const double D[2][4] = {{0.0164, 0.0085, 0.00832, 0.00821},
                        {0.0235, 0.0121, 0.0119, 0.0116}};
const double Sa[2][4] = {{0.00139, 0.000218, 0.00197, 0.0106},
                         {0.00977, 0.162, 0.156, 0.535}};
const double Sr[2][4] = {{1.77139, 0.533218, 3.31197, 0.0106},
                         {1.23977, 0.529, 2.436, 0.535}};
const double Sf[2][4] = {{0.0, 0.0, 0.0, 0.0}, {0.00395, 0.0262, 0.0718, 0.346}};
const double nu[2][4] = {{0.0, 0.0, 0.0, 0.0}, {2.49, 2.43, 2.42, 2.42}};
const double chi[2][4] = {{0.0, 0.0, 0.0, 0.0}, {0.9675, 0.03250, 0.0, 0.0}};
const double Ss[2][4][4] = {{{ 0.0,   0.0,  0.0, 0.0},
                             {1.77,   0.0,  0.0, 0.0},
                             { 0.0, 0.533,  0.0, 0.0},
                             { 0.0,   0.0, 3.31, 0.0}},
                            {{ 0.0,   0.0,  0.0, 0.0},
                             {1.23,   0.0,  0.0, 0.0},
                             { 0.0, 0.367,  0.0, 0.0},
                             { 0.0,   0.0, 2.28, 0.0}}};

// initial eigenvalue approximation
double k_eff = 1.0;

//////  Bilinear and linear forms - axisymmetric arrangement  ////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar int_x_u_v(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (e->x[i] * u->val[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_x_grad_u_grad_v(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * e->x[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

//////////   Eq 1   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar biform_0_0(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (D[e->marker - 1][0]) * int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) +
         (Sr[e->marker - 1][0]) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_surf_0_0(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (0.5) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar liform_0(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (chi[e->marker - 1][0] / k_eff) * (nu[e->marker - 1][0] * Sf[e->marker - 1][0] * ext->fn[0]->val[i] +
                                                   nu[e->marker - 1][1] * Sf[e->marker - 1][1] * ext->fn[1]->val[i] +
                                                   nu[e->marker - 1][2] * Sf[e->marker - 1][2] * ext->fn[2]->val[i] +
                                                   nu[e->marker - 1][3] * Sf[e->marker - 1][3] * ext->fn[3]->val[i])
                                      * e->x[i] * v->val[i];
  return result;
}

//////////   Eq 2   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar biform_1_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (D[e->marker - 1][1]) * int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) +
         (Sr[e->marker - 1][1]) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_surf_1_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (0.5) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_1_0(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (- Ss[e->marker - 1][1][0]) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar liform_1(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (chi[e->marker - 1][1] / k_eff) * (nu[e->marker - 1][0] * Sf[e->marker - 1][0] * ext->fn[0]->val[i] +
                                                   nu[e->marker - 1][1] * Sf[e->marker - 1][1] * ext->fn[1]->val[i] +
                                                   nu[e->marker - 1][2] * Sf[e->marker - 1][2] * ext->fn[2]->val[i] +
                                                   nu[e->marker - 1][3] * Sf[e->marker - 1][3] * ext->fn[3]->val[i])
                                      * e->x[i] * v->val[i];
  return result;
}

//////////   Eq 3   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar biform_2_2(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (D[e->marker - 1][2]) * int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) +
         (Sr[e->marker - 1][2]) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_surf_2_2(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (0.5) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_2_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (- Ss[e->marker - 1][2][1]) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar liform_2(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (chi[e->marker - 1][2] / k_eff) * (nu[e->marker - 1][0] * Sf[e->marker - 1][0] * ext->fn[0]->val[i] +
                                                   nu[e->marker - 1][1] * Sf[e->marker - 1][1] * ext->fn[1]->val[i] +
                                                   nu[e->marker - 1][2] * Sf[e->marker - 1][2] * ext->fn[2]->val[i] +
                                                   nu[e->marker - 1][3] * Sf[e->marker - 1][3] * ext->fn[3]->val[i])
                                      * e->x[i] * v->val[i];
  return result;
}

//////////   Eq 4   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar biform_3_3(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (D[e->marker - 1][3]) * int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) +
         (Sr[e->marker - 1][3]) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_surf_3_3(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (0.5) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_3_2(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (- Ss[e->marker - 1][3][2]) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar liform_3(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (chi[e->marker - 1][3] / k_eff) * (nu[e->marker - 1][0] * Sf[e->marker - 1][0] * ext->fn[0]->val[i] +
                                                   nu[e->marker - 1][1] * Sf[e->marker - 1][1] * ext->fn[1]->val[i] +
                                                   nu[e->marker - 1][2] * Sf[e->marker - 1][2] * ext->fn[2]->val[i] +
                                                   nu[e->marker - 1][3] * Sf[e->marker - 1][3] * ext->fn[3]->val[i])
                                      * e->x[i] * v->val[i];
  return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void source_fn(int n, scalar* a, scalar* b, scalar* c, scalar* d, scalar* out)
{
  for (int i = 0; i < n; i++)
  {
    out[i] = (nu[1][0] * Sf[1][0] * a[i] +
        nu[1][1] * Sf[1][1] * b[i] +
        nu[1][2] * Sf[1][2] * c[i] +
        nu[1][3] * Sf[1][3] * d[i]);
  }
}

double integrate(MeshFunction* sln, int marker)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);

  double integral = 0.0;
  Element* e;
  Mesh* mesh = sln->get_mesh();

  for_all_active_elements(e, mesh)
  {
    if (e->marker == marker)
    {
      update_limit_table(e->get_mode());
      sln->set_active_element(e);
      RefMap* ru = sln->get_refmap();
      int o = sln->get_fn_order() + ru->get_inv_ref_order();
      limit_order(o);
      sln->set_quad_order(o, FN_VAL);
      scalar *uval = sln->get_fn_values();
      double* x = ru->get_phys_x(o);
      double result = 0.0;
      h1_integrate_expression(x[i] * uval[i]);
      integral += result;
    }
  }

  return 2.0 * M_PI * integral;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{

  Mesh mesh;
  H2DReader mloader;
  mloader.load("reactor.mesh", &mesh);
  for (int i = 0; i < init_refinements; i++) mesh.refine_all_elements();

  H1Shapeset shapeset;
  PrecalcShapeset pss1(&shapeset);
  PrecalcShapeset pss2(&shapeset);
  PrecalcShapeset pss3(&shapeset);
  PrecalcShapeset pss4(&shapeset);

  Solution sln1, sln2, sln3, sln4;
  Solution iter1, iter2, iter3, iter4;
  Solution ref1, ref2, ref3, ref4;
  iter1.set_const(&mesh, 1.00);
  iter2.set_const(&mesh, 1.00);
  iter3.set_const(&mesh, 1.00);
  iter4.set_const(&mesh, 1.00);
  UmfpackSolver umfpack;

  H1Space space1(&mesh, &shapeset);
  H1Space space2(&mesh, &shapeset);
  H1Space space3(&mesh, &shapeset);
  H1Space space4(&mesh, &shapeset);
  space1.set_bc_types(bc_types);
  space2.set_bc_types(bc_types);
  space3.set_bc_types(bc_types);
  space4.set_bc_types(bc_types);
  space1.set_uniform_order(init_order);
  space2.set_uniform_order(init_order);
  space3.set_uniform_order(init_order);
  space4.set_uniform_order(init_order);

  WeakForm wf(4);
  wf.add_biform(0, 0, callback(biform_0_0));
  wf.add_biform(1, 1, callback(biform_1_1));
  wf.add_biform(1, 0, callback(biform_1_0));
  wf.add_biform(2, 2, callback(biform_2_2));
  wf.add_biform(2, 1, callback(biform_2_1));
  wf.add_biform(3, 3, callback(biform_3_3));
  wf.add_biform(3, 2, callback(biform_3_2));

  wf.add_liform(0, callback(liform_0), marker_core, 4, &iter1, &iter2, &iter3, &iter4);
  wf.add_liform(1, callback(liform_1), marker_core, 4, &iter1, &iter2, &iter3, &iter4);
  wf.add_liform(2, callback(liform_2), marker_core, 4, &iter1, &iter2, &iter3, &iter4);
  wf.add_liform(3, callback(liform_3), marker_core, 4, &iter1, &iter2, &iter3, &iter4);

  wf.add_biform_surf(0, 0, callback(biform_surf_0_0), bc_vacuum);
  wf.add_biform_surf(1, 1, callback(biform_surf_1_1), bc_vacuum);
  wf.add_biform_surf(2, 2, callback(biform_surf_2_2), bc_vacuum);
  wf.add_biform_surf(3, 3, callback(biform_surf_3_3), bc_vacuum);

  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(4, &space1, &space2, &space3, &space4);
  sys.set_pss(4, &pss1, &pss2, &pss3, &pss4);

  ScalarView view1("Neutron flux 1", 0, 0, 320, 600);
  ScalarView view2("Neutron flux 2", 350, 0, 320, 600);
  ScalarView view3("Neutron flux 3", 700, 0, 320, 600);
  ScalarView view4("Neutron flux 4", 1050, 0, 320, 600);
  view1.show_mesh(false);
  view2.show_mesh(false);
  view3.show_mesh(false);
  view4.show_mesh(false);

  int ndofs = 0;
  ndofs += space1.assign_dofs(ndofs);
  ndofs += space2.assign_dofs(ndofs);
  ndofs += space3.assign_dofs(ndofs);
  ndofs += space4.assign_dofs(ndofs);

  // Main power iterations loop
  bool eigen_done = false; int it = 0;
  do
  {
    info("\n------------ Eigen Iteration %d ----------------------\n", it++);

    sys.assemble();
    sys.solve(4, &sln1, &sln2, &sln3, &sln4);

    // visualization
    view1.show(&sln1);    view2.show(&sln2);
    view3.show(&sln3);    view4.show(&sln4);


    SimpleFilter source(source_fn, &sln1, &sln2, &sln3, &sln4);
    SimpleFilter source_prev(source_fn, &iter1, &iter2, &iter3, &iter4);

    // compute eigenvalue
    double k_new = k_eff * (integrate(&source, marker_core) / integrate(&source_prev, marker_core));
    info("Greatest eigenvalue approximation: %g, (relative error: %g)", k_new, fabs((k_eff - k_new) / k_new));

    // stopping criterion
    if (fabs((k_eff - k_new) / k_new) < 1e-5) eigen_done = true;

    // update eigenvectors
    iter1.copy(&sln1);    iter2.copy(&sln2);
    iter3.copy(&sln3);    iter4.copy(&sln4);

    // update eigenvalue
    k_eff = k_new;
  }
  while (!eigen_done);

  View::wait();
  return 0;
}
