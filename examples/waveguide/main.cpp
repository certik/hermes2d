#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This example solves adaptively the electric field in a simplified microwave oven.
// The waves are generated using a harmonic surface current on the right-most edge.
// (Such small cavity is present in every microwave oven). There is a circular
// load located in the middle of the main cavity, defined through a different
// permittivity -- see function in_load(...). One can either use a mesh that is
// aligned to the load via curvilinear elements (ALIGN_MESH = true), or an unaligned
// mesh (ALIGN_MESH = false). Convergence graphs are saved both wrt. the dof number
// and cpu time.
//
// PDE: time-harmonic Maxwell's equations;
//      there is circular load in the middle of the large cavity, whose permittivity
//      is different from the rest of the domain.
//
// Domain: square cavity with another small square cavity attached from outside
//         on the right.
//
// Meshes: you can either use "oven_load_circle.mesh" containing curved elements
//         aligned with the circular load, or "oven_load_square.mesh" which is not
//         aligned.
//
// BC: perfect conductor on the boundary except for the right-most edge of the small
//     cavity, where a harmonic surface current is prescribed
//
// The following parameters can be changed:

const bool SOLVE_ON_COARSE_MESH = false; // If true, coarse mesh FE problem is solved in every adaptivity step.
                                         // If false, projection of the fine mesh solution on the coarse mesh is used. 
const int INIT_REF_NUM = 0;              // Number of initial uniform mesh refinements.
const int P_INIT = 2;                    // Initial polynomial degree. NOTE: The meaning is different from
                                         // standard continuous elements in the space H1. Here, P_INIT refers
                                         // to the maximum poly order of the tangential component, and polynomials
                                         // of degree P_INIT + 1 are present in element interiors. P_INIT = 0
                                         // is for Whitney elements.
const bool ALIGN_MESH = true;           // if ALIGN_MESH == true, curvilinear elements aligned with the
                                         // circular load are used, otherwise one uses a non-aligned mesh.
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See the User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 2.0;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows
                                         // over this limit. This is to prevent h-adaptivity to go on forever.

// Problem parameters.
const double e_0   = 8.8541878176 * 1e-12;
const double mu_0   = 1.256 * 1e-6;
const double e_r = 1.0;
const double mu_r = 1.0;
const double rho = 3820.0;
const double Cp = 7.531000;
const double freq = 1.0*2450000000.0;
const double omega = 2 * M_PI * freq;
const double c = 1 / sqrt(e_0 * mu_0);
const double kappa  = 2 * M_PI * freq * sqrt(e_0 * mu_0);
const double J = 0.0000033333;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == 2) return BC_ESSENTIAL; // perfect conductor BC
  else return BC_NATURAL;               // impedance BC
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

// Geometry of the load.
bool in_load(double x, double y)
{
  double cx = -0.152994121;
  double cy =  0.030598824;
  double r = 0.043273273;
  if (sqr(cx - x) + sqr(cy - y) < sqr(r)) return true;
  else return false;
}

// Gamma as a function of x, y.
double gam(int marker, double x, double y)
{
  if (ALIGN_MESH && marker == 1) return 0.03;
  if (!ALIGN_MESH && in_load(x,y)) {
    double cx = -0.152994121;  double cy =  0.030598824;
    double r = sqrt(sqr(cx - x) + sqr(cy - y));
    return (0.03 + 1)/2.0 - (0.03 - 1) * atan(10.0*(r -  0.043273273)) / M_PI;
  }
  return 0.0;
}
double gam(int marker, Ord x, Ord y)
{  return 0.0; }


// Relative permittivity as a function of x, y.
double er(int marker, double x, double y)
{
  if (ALIGN_MESH && marker == 1) return 7.5;
  if (!ALIGN_MESH && in_load(x,y)) {
    double cx = -0.152994121;  double cy =  0.030598824;
    double r = sqrt(sqr(cx - x) + sqr(cy - y));
    return (7.5 + 1)/2.0 - (7.5 - 1) * atan(10.0*(r -  0.043273273)) / M_PI;
  }
  return 1.0;
}
double er(int marker, Ord x, Ord y)
{  return 1.0; }

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  if (ALIGN_MESH) mloader.load("oven_load_circle.mesh", &mesh);
  else mloader.load("oven_load_square.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Create an Hcurl space.
  HcurlSpace space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form_surf(callback(linear_form_surf));

  // Initialize views.
  VectorView eview("Electric field", 0, 0, 580, 400);
  OrderView ord("Order", 590, 0, 550, 400);

  /*
  // View the basis functions.
  VectorBaseView bview;
  vbview.show(&space);
  */

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Initialize refinements selector.
  HcurlProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize the coarse mesh problem.
  LinSystem ls(&wf, &space);

  // Adaptivity loop:
  int as = 1; bool done = false;
  Solution sln_coarse, sln_fine;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Assemble and solve the fine mesh problem.
    info("Solving on fine mesh.");
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(&sln_fine);

    // Either solve on coarse mesh or project the fine mesh solution 
    // on the coarse mesh.
    if (SOLVE_ON_COARSE_MESH) {
      info("Solving on coarse mesh.");
      ls.assemble();
      ls.solve(&sln_coarse);
    }
    else {
      info("Projecting fine mesh solution on coarse mesh.");
      int proj_type = 2;    // Hcurl projection.
      ls.project_global(&sln_fine, &sln_coarse, proj_type);
    }

    // Time measurement.
    cpu_time.tick();

    // Show real part of the solution.
    AbsFilter abs(&sln_coarse);
    eview.set_min_max_range(0, 4e3);
    eview.show(&abs);
    ord.show(&space);

    // Skip visualization time. 
    cpu_time.tick(H2D_SKIP);

    // Calculate error estimate wrt. fine mesh solution.
    info("Calculating error.");
    HcurlAdapt hp(&ls);
    hp.set_solutions(&sln_coarse, &sln_fine);
    hp.set_biform(callback(hcurl_form_kappa));  // Adaptivity is done using an "energetic norm".
    double err_est_adapt = hp.calc_error() * 100;
    double err_est_hcurl = hcurl_error(&sln_coarse, &sln_fine) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est: %g%%", 
         ls.get_num_dofs(), rs.get_num_dofs(), err_est_hcurl);

    // Add entries to DOF convergence graphs.
    graph_dof_est.add_values(ls.get_num_dofs(), err_est_hcurl);
    graph_dof_est.save("conv_dof_est.dat");

    // Add entries to CPU convergence graphs.
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_hcurl);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est_adapt too large, adapt the mesh.
    if (err_est_hcurl < ERR_STOP) done = true;
    else {
      info("Adapting coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      if (ls.get_num_dofs() >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the fine mesh solution - the final result.
  eview.set_title("Final solution");
  eview.show(&sln_fine);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

