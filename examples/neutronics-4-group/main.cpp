#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"

#include "hermes2d.h"

// This example solves a 4-group neutron diffusion equation in the reactor core.
// The eigenproblem is solved using power interations.
//
// The reactor neutronics is given by the following eigenproblem:
//
//  - \nabla \cdot D_g \nabla \phi_g + \Sigma_{Rg}\phi_g - \sum_{g' \neq g} \Sigma_s^{g'\to g} \phi_{g'} =
//  = \frac{\chi_g}{k_{eff}} \sum_{g'} \nu_{g'} \Sigma_{fg'}\phi_{g'}
//
// where 1/k_{eff} is eigenvalue and \phi_g, g = 1,...,4 are eigenvectors (neutron fluxes). The current problem
// is posed in a 3D cylindrical axisymmetric geometry, leading to a 2D problem with r-z as the independent spatial 
// coordinates. The corresponding diffusion operator is given by (r = x, z = y):
//
//	\nabla \cdot D \nabla \phi = \frac{1}{x} (x D \phi_x)_x  + (D \phi_y)_y 
//
// BC:
//
// Homogeneous neumann on symmetry axis,
// d \phi_g / d n = - 0.5 \phi_g   elsewhere
//
// The eigenproblem is numerically solved using common technique known as the power method (power iterations):
//
//  1) Make an initial estimate of \phi_g and k_{eff}
//  2) For n = 1, 2,...
//         solve for \phi_g using previous k_prev
//         solve for new k_{eff}
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{new}
//               k_new =  k_prev -------------------------------------------------------------------------
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{prev}
//  3) Stop iterations when
//
//     |   k_new - k_prev  |
//     | ----------------- |  < epsilon
//     |       k_new       |
//
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 2;                  // Number of initial uniform mesh refinements.
const int P_INIT = 2;                        // Initial polynomial degree of all mesh elements.
const double ERROR_STOP = 1e-5;              // Tolerance for the eigenvalue.

// Element markers.
const int MAT_REFLECTOR = 1;
const int MAT_CORE = 2;

// Boundary merkers.
const int BDY_VACUUM = 1;
const int BDY_SYMMETRY = 2;

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

// Reflector properties (0) core properties (1).
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

// Initial eigenvalue approximation.
double k_eff = 1.0;

// Weak forms.
#include "forms.cpp"

// Source function.
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

// Integral over the active core.
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

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("reactor.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Solution variables.
  Solution sln1, sln2, sln3, sln4;
  Solution iter1, iter2, iter3, iter4;

  // Define initial conditions.
  info("Setting initial conditions.");
  iter1.set_const(&mesh, 1.00);
  iter2.set_const(&mesh, 1.00);
  iter3.set_const(&mesh, 1.00);
  iter4.set_const(&mesh, 1.00);

  // Create H1 spaces with default shapesets.
  H1Space space1(&mesh, bc_types, essential_bc_values, P_INIT);
  H1Space space2(&mesh, bc_types, essential_bc_values, P_INIT); 
  H1Space space3(&mesh, bc_types, essential_bc_values, P_INIT); 
  H1Space space4(&mesh, bc_types, essential_bc_values, P_INIT); 

  // Initialize the weak formulation.
  WeakForm wf(4);
  wf.add_matrix_form(0, 0, callback(biform_0_0));
  wf.add_matrix_form(1, 1, callback(biform_1_1));
  wf.add_matrix_form(1, 0, callback(biform_1_0));
  wf.add_matrix_form(2, 2, callback(biform_2_2));
  wf.add_matrix_form(2, 1, callback(biform_2_1));
  wf.add_matrix_form(3, 3, callback(biform_3_3));
  wf.add_matrix_form(3, 2, callback(biform_3_2));
  wf.add_vector_form(0, callback(liform_0), MAT_CORE, Tuple<MeshFunction*>(&iter1, &iter2, &iter3, &iter4));
  wf.add_vector_form(1, callback(liform_1), MAT_CORE, Tuple<MeshFunction*>(&iter1, &iter2, &iter3, &iter4));
  wf.add_vector_form(2, callback(liform_2), MAT_CORE, Tuple<MeshFunction*>(&iter1, &iter2, &iter3, &iter4));
  wf.add_vector_form(3, callback(liform_3), MAT_CORE, Tuple<MeshFunction*>(&iter1, &iter2, &iter3, &iter4));
  wf.add_matrix_form_surf(0, 0, callback(biform_surf_0_0), BDY_VACUUM);
  wf.add_matrix_form_surf(1, 1, callback(biform_surf_1_1), BDY_VACUUM);
  wf.add_matrix_form_surf(2, 2, callback(biform_surf_2_2), BDY_VACUUM);
  wf.add_matrix_form_surf(3, 3, callback(biform_surf_3_3), BDY_VACUUM);

  // Initialize coarse mesh problem.
  LinSystem ls(&wf, Tuple<Space*>(&space1, &space2, &space3, &space4));

  // Initialize views.
  ScalarView view1("Neutron flux 1", 0, 0, 320, 600);
  ScalarView view2("Neutron flux 2", 350, 0, 320, 600);
  ScalarView view3("Neutron flux 3", 700, 0, 320, 600);
  ScalarView view4("Neutron flux 4", 1050, 0, 320, 600);

  // Do not show meshes.
  view1.show_mesh(false);
  view2.show_mesh(false);
  view3.show_mesh(false);
  view4.show_mesh(false);

  // Main power iteration loop:
  int iter = 0; bool done = false; 
  bool rhs_only = false;
  do
  {
    info("------------ Power iteration %d:", iter);

    // Assemble and solve the coarse mesh problem.
    ls.assemble(rhs_only);
    ls.solve(Tuple<Solution*>(&sln1, &sln2, &sln3, &sln4));

    // Show solutions.
    view1.show(&sln1);    
    view2.show(&sln2);
    view3.show(&sln3);    
    view4.show(&sln4);

    SimpleFilter source(source_fn, &sln1, &sln2, &sln3, &sln4);
    SimpleFilter source_prev(source_fn, &iter1, &iter2, &iter3, &iter4);

    // Compute eigenvalue.
    double k_new = k_eff * (integrate(&source, MAT_CORE) / integrate(&source_prev, MAT_CORE));
    info("ndof: %d, %d, %d, %d", space1.get_num_dofs(),space2.get_num_dofs(), 
                                  space3.get_num_dofs(), space4.get_num_dofs());
    info("Largest eigenvalue: (%.8g, %.8g), rel error: %g", k_eff, k_new, fabs((k_eff - k_new) / k_new));

    // Stopping criterion.
    if (fabs((k_eff - k_new) / k_new) < ERROR_STOP) done = true;

    // Save solutions for the next iteration.
    iter1.copy(&sln1);    
    iter2.copy(&sln2);
    iter3.copy(&sln3);    
    iter4.copy(&sln4);

    // Update eigenvalue.
    k_eff = k_new;
    
    // Don't need to reassemble the system matrix in further iterations,
    // only the rhs changes to reflect the progressively updated source.
    rhs_only = true;

    iter++;
  }
  while (!done);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
