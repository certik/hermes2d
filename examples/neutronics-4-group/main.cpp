#include "hermes2d.h"
#include "solver_umfpack.h"

// This example solves a 4-group neutron diffusion equation in the reactor core.
// The eigenproblem is solved using power interations.
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


const int P_INIT = 1;
const int INIT_REF_NUM = 4;

// Area markers
const int marker_reflector = 1;
const int marker_core = 2;

// Boundary indices
const int bc_vacuum = 1;
const int bc_sym = 2;

// Boundary condition types
int bc_types(int marker)
{
  return BC_NATURAL;
}

// Reflector properties (0) core properties (1),
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

// Initial eigenvalue approximation
double k_eff = 1.0;

// Weak forms
#include "forms.cpp"

// source function
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

// Integral over the active core
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
      sln->set_quad_order(o, H2D_FN_VAL);
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
  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("reactor.mesh", &mesh);

  // initial uniform refinements
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss1(&shapeset);
  PrecalcShapeset pss2(&shapeset);
  PrecalcShapeset pss3(&shapeset);
  PrecalcShapeset pss4(&shapeset);

  // solution variables
  Solution sln1, sln2, sln3, sln4;
  Solution iter1, iter2, iter3, iter4;
  Solution ref1, ref2, ref3, ref4;
  iter1.set_const(&mesh, 1.00);
  iter2.set_const(&mesh, 1.00);
  iter3.set_const(&mesh, 1.00);
  iter4.set_const(&mesh, 1.00);

  // matrix solver
  UmfpackSolver umfpack;

  // create finite element spaces
  H1Space space1(&mesh, &shapeset);
  H1Space space2(&mesh, &shapeset);
  H1Space space3(&mesh, &shapeset);
  H1Space space4(&mesh, &shapeset);
  space1.set_bc_types(bc_types);
  space2.set_bc_types(bc_types);
  space3.set_bc_types(bc_types);
  space4.set_bc_types(bc_types);
  space1.set_uniform_order(P_INIT);
  space2.set_uniform_order(P_INIT);
  space3.set_uniform_order(P_INIT);
  space4.set_uniform_order(P_INIT);

  // initialize the weak formulation
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

  // initialize the LinSystem class
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(4, &space1, &space2, &space3, &space4);
  sys.set_pss(4, &pss1, &pss2, &pss3, &pss4);

  // visualization
  ScalarView view1("Neutron flux 1", 0, 0, 320, 600);
  ScalarView view2("Neutron flux 2", 350, 0, 320, 600);
  ScalarView view3("Neutron flux 3", 700, 0, 320, 600);
  ScalarView view4("Neutron flux 4", 1050, 0, 320, 600);
  view1.show_mesh(false);
  view2.show_mesh(false);
  view3.show_mesh(false);
  view4.show_mesh(false);

  // enumerate basis functions
  int ndof = assign_dofs(4, &space1, &space2, &space3, &space4);

  // Main power iteration loop
  bool eigen_done = false; int it = 0;
  do
  {
    info("\n------------ Power iteration %d ----------------------\n", it++);

    sys.assemble();
    sys.solve(4, &sln1, &sln2, &sln3, &sln4);

    // visualization
    view1.show(&sln1);    view2.show(&sln2);
    view3.show(&sln3);    view4.show(&sln4);

    SimpleFilter source(source_fn, &sln1, &sln2, &sln3, &sln4);
    SimpleFilter source_prev(source_fn, &iter1, &iter2, &iter3, &iter4);

    // compute eigenvalue
    double k_new = k_eff * (integrate(&source, marker_core) / integrate(&source_prev, marker_core));
    info("Largest eigenvalue (est): %g, rel error: %g", k_new, fabs((k_eff - k_new) / k_new));

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
