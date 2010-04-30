#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

//  This example has a known exact solution. It describes an electromagnetic wave that hits
//  a screen under the angle of 45 degrees, causing a singularity at the tip of the screen.
//  Convergence graphs saved (both exact error and error estimate, and both wrt. dof number
//  and cpu time).
//
//  PDE: time-harmonic Maxwell's equations
//
//  Known exact solution, see the function exact()
//
//  Domain: square domain cut from the midpoint of the left edge to the center (center
//          point of left edge duplicated)
//
//  Meshes: you can either use "screen-quad.mesh" (quadrilateral mesh) or
//          "screen-tri.mesh" (triangular mesh). See the command mesh.load(...) below
//
//  BC: tangential component of solution taken from known exact solution (essential BC),
//      see function bc_values(...) below
//
// The following parameters can be changed:

const int P_INIT = 0;             // Initial polynomial degree. NOTE: The meaning is different from
                                  // standard continuous elements in the space H1. Here, P_INIT refers
                                  // to the maximum poly order of the tangential component, and polynomials
                                  // of degree P_INIT + 1 are present in element interiors. P_INIT = 0
                                  // is for Whitney elements.
const double THRESHOLD = 0.5;     // This is a quantitative parameter of the adapt(...) function and
                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;           // Adaptive strategy:
                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                  //   error is processed. If more elements have similar errors, refine
                                  //   all to keep the mesh symmetric.
                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                  //   than THRESHOLD times maximum element error.
                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                  //   than THRESHOLD.
                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int ADAPT_TYPE = 0;         // Type of automatic adaptivity:
                                  // ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                  // ADAPT_TYPE = 1 ... adaptive h-FEM,
                                  // ADAPT_TYPE = 2 ... adaptive p-FEM.
const bool ISO_ONLY = false;      // Isotropic refinement flag (concerns quadrilateral elements only).
                                  // ISO_ONLY = false ... anisotropic refinement of quad elements
                                  // is allowed (default),
                                  // ISO_ONLY = true ... only isotropic refinements of quad elements
                                  // are allowed.
const int MESH_REGULARITY = -1;   // Maximum allowed level of hanging nodes:
                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                  // Note that regular meshes are not supported, this is due to
                                  // their notoriously bad performance.
const double ERR_STOP = 0.1;      // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 50000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// problem constants
const double e_0  = 8.8541878176 * 1e-12;
const double mu_0 = 1.256 * 1e-6;
const double k = 1.0;

// exact solution
#include "exact_sol.cpp"

// boundary conditions
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

double2 tau[5] = { { 0, 0}, { 1, 0 },  { 0, 1 }, { -1, 0 }, { 0, -1 } };

cplx bc_values(int marker, double x, double y)
{
  scalar dx, dy;
  return exact0(x, y, dx, dy)*tau[marker][0] + exact1(x, y, dx, dy)*tau[marker][1];
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_curl_e_curl_f<Real, Scalar>(n, wt, u, v) - int_e_f<Real, Scalar>(n, wt, u, v);
}

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("screen-quad.mesh", &mesh);
//    mloader.load("screen-tri.mesh", &mesh);

  // initialize the shapeset and the cache
  HcurlShapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create finite element space
  HcurlSpace space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // enumerate degrees of freedom
  int ndof = assign_dofs(&space);

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form), H2D_SYM);

  // visualize solution and mesh
  ScalarView Xview_r("Electric field X - real",   0, 0, 450, 420);
  ScalarView Yview_r("Electric field Y - real", 460, 0, 450, 420);
  ScalarView Xview_i("Electric field X - imag", 920, 0, 450, 420);
  ScalarView Yview_i("Electric field Y - imag", 1380, 0, 450, 420);
  OrderView  ord("Polynomial Orders", 0, 460, 450, 420);

  /*
  // view the basis functions
  VectorBaseView bview;
  vbview.show(&space);
  vbview.wait_for_keypress();
  */

  // matrix solver
  UmfpackSolver solver;

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof_est, graph_dof_exact, graph_cpu_est, graph_cpu_exact;

  // adaptivity loop
  int it = 1;
  bool done = false;
  TimePeriod cpu_time;
  Solution sln_coarse, sln_fine;
  do
  {
    info("!---- Adaptivity step %d ---------------------------------------------", it); it++;

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // solve the coarse mesh problem
    LinSystem sys(&wf, &solver);
    sys.set_spaces(1, &space);
    sys.set_pss(1, &pss);
    sys.assemble();
    sys.solve(1, &sln_coarse);

    // time measurement
    cpu_time.tick();

    // calculating error wrt. exact solution
    Solution ex;
    ex.set_exact(&mesh, exact);
    double err_exact = 100 * hcurl_error(&sln_coarse, &ex);
    info("Exact solution error: %g%%", err_exact);

    // visualization
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

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // solve the fine mesh problem
    RefSystem ref(&sys);
    ref.assemble();
    ref.solve(1, &sln_fine);

    // calculate error estimate wrt. fine mesh solution
    HcurlOrthoHP hp(1, &space);
    double err_est_adapt = hp.calc_error(&sln_coarse, &sln_fine) * 100;
    double err_est_hcurl = hcurl_error(&sln_coarse, &sln_fine) * 100;
    
    // report results
    cpu_time.tick();
    info("Error estimate (adapt): %g%%", err_est_adapt);
    info("Error estimate (hcurl): %g%%", err_est_hcurl);

    // add entries to DOF convergence graphs
    graph_dof_exact.add_values(space.get_num_dofs(), err_exact);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_dof_est.add_values(space.get_num_dofs(), err_est_hcurl);
    graph_dof_est.save("conv_dof_est.dat");

    // add entries to CPU convergence graphs
    graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact);
    graph_cpu_exact.save("conv_cpu_exact.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_hcurl);
    graph_cpu_est.save("conv_cpu_est.dat");

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // if err_est_adapt too large, adapt the mesh
    if (err_est_adapt < ERR_STOP) done = true;
    else {
      hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE, ISO_ONLY, MESH_REGULARITY);
      ndof = assign_dofs(&space);
      if (ndof >= NDOF_STOP) done = true;
    }

    // time measurement
    cpu_time.tick();
  }
  while (!done);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // wait for all views to be closed
  View::wait();
  return 0;
}

