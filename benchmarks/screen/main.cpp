#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This benchmark problem with known exact solution describes an electromagnetic wave that hits
//  a screen under the angle of 45 degrees, causing a very strong singularity at the tip of the 
//  screen. Convergence graphs saved (both exact error and error estimate, and both wrt. dof number.
//  and cpu time).
//
//  PDE: time-harmonic Maxwell's equations.
//
//  Known exact solution, see the function exact().
//
//  Domain: square domain cut from the midpoint of the left edge to the center (center
//          point of left edge duplicated).
//
//  Meshes: you can either use "screen-quad.mesh" (quadrilateral mesh) or
//          "screen-tri.mesh" (triangular mesh). See the command mesh.load(...) below.
//
//  BC: tangential component of solution taken from known exact solution (essential BC),
//      see function essential_bc_values() below.
//
// The following parameters can be changed:

const bool SOLVE_ON_COARSE_MESH = false; // If true, coarse mesh FE problem is solved in every adaptivity step.
                                         // If false, projection of the fine mesh solution on the coarse mesh is used. 
const int INIT_REF_NUM = 1;              // Number of initial uniform mesh refinements.
const int P_INIT = 1;                    // Initial polynomial degree. NOTE: The meaning is different from
                                         // standard continuous elements in the space H1. Here, P_INIT refers
                                         // to the maximum poly order of the tangential component, and polynomials
                                         // of degree P_INIT + 1 are present in element interiors. P_INIT = 0
                                         // is for Whitney elements.
const double THRESHOLD = 0.5;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO_H; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
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
const int NDOF_STOP = 50000;             // Adaptivity process stops when the number of degrees of freedom grows
                                         // over this limit. This is to prevent h-adaptivity to go on forever.

// Problem parameters.
const double e_0  = 8.8541878176 * 1e-12;
const double mu_0 = 1.256 * 1e-6;
const double k = 1.0;

// Exact solution.
#include "exact_solution.cpp"

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Unit tangential vectors to the boundary. 
double2 tau[5] = { { 0, 0}, { 1, 0 },  { 0, 1 }, { -1, 0 }, { 0, -1 } };

// Essential boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  scalar dx, dy;
  return exact0(x, y, dx, dy)*tau[ess_bdy_marker][0] + exact1(x, y, dx, dy)*tau[ess_bdy_marker][1];
}

// Weak forms.
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_curl_e_curl_f<Real, Scalar>(n, wt, u, v) - int_e_f<Real, Scalar>(n, wt, u, v);
}

int main(int argc, char* argv[])
{
  // Time measurement
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("screen-quad.mesh", &mesh);    // quadrilaterals
  // mloader.load("screen-tri.mesh", &mesh);  // triangles

  // Perform initial mesh refinements.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Create an Hcurl space with default shapeset.
  HcurlSpace space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form), H2D_SYM);

  // Initialize views.
  ScalarView Xview_r("Electric field X - real",   0, 0, 300, 280);
  ScalarView Yview_r("Electric field Y - real", 310, 0, 300, 280);
  ScalarView Xview_i("Electric field X - imag", 620, 0, 300, 280);
  ScalarView Yview_i("Electric field Y - imag", 930, 0, 300, 280);
  OrderView  ord("Polynomial Orders", 0, 335, 300, 280);

  /*
  // View the basis functions.
  VectorBaseView bview;
  vbview.show(&space);
  View::wait(H2DV_WAIT_KEYPRESS);
  */

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_dof_exact, graph_cpu_est, graph_cpu_exact;

  // Initialize refinement selector.
  HcurlProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize the coarse mesh problem.
  LinSystem ls(&wf, &space);

  // Adaptivity loop.
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

    // Calculate error wrt. exact solution.
    info("Calculating error (exact).");
    Solution ex;
    ex.set_exact(&mesh, exact);
    double err_exact = 100 * hcurl_error(&sln_coarse, &ex);

    // Visualization.
    RealFilter real(&sln_coarse);
    ImagFilter imag(&sln_coarse);
    Xview_r.set_min_max_range(-3.0, 1.0);
    //Xview_r.show_scale(false);
    Xview_r.show(&real, H2D_EPS_NORMAL, H2D_FN_VAL_0);
    Yview_r.set_min_max_range(-4.0, 4.0);
    //Yview_r.show_scale(false);
    Yview_r.show(&real, H2D_EPS_NORMAL, H2D_FN_VAL_1);
    Xview_i.set_min_max_range(-1.0, 4.0);
    //Xview_i.show_scale(false);
    Xview_i.show(&imag, H2D_EPS_NORMAL, H2D_FN_VAL_0);
    Yview_i.set_min_max_range(-4.0, 4.0);
    //Yview_i.show_scale(false);
    Yview_i.show(&imag, H2D_EPS_NORMAL, H2D_FN_VAL_1);
    ord.show(&space);

    // Skip exact error calculation and visualization time. 
    cpu_time.tick(H2D_SKIP);

    // Calculate error estimate wrt. fine mesh solution.
    info("Calculating error (est).");
    HcurlAdapt hp(&ls);
    hp.set_solutions(&sln_coarse, &sln_fine);
    double err_est_adapt = hp.calc_error() * 100;
    double err_est_hcurl = hcurl_error(&sln_coarse, &sln_fine) * 100;
    
    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est: %g%%, err_exact: %g%%", 
         ls.get_num_dofs(), rs.get_num_dofs(), err_est_hcurl, err_exact);

    // Add entries to DOF convergence graphs.
    graph_dof_exact.add_values(ls.get_num_dofs(), err_exact);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_dof_est.add_values(ls.get_num_dofs(), err_est_hcurl);
    graph_dof_est.save("conv_dof_est.dat");

    // Add entries to CPU convergence graphs.
    graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact);
    graph_cpu_exact.save("conv_cpu_exact.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_hcurl);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est_adapt too large, adapt the mesh.
    if (err_est_adapt < ERR_STOP) done = true;
    else {
      info("Adapting coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      if (ls.get_num_dofs() >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (!done);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

