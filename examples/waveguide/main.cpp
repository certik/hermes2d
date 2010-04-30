#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

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
//      is different from the rest of the domain
//
// Domain: square cavity with another small square cavity attached from outside
//         on the right
//
// Meshes: you can either use "oven_load_circle.mesh" containing curved elements
//         aligned with the circular load, or "oven_load_square.mesh" which is not
//         aligned.
//
// BC: perfect conductor on the boundary except for the right-most edge of the small
//     cavity, where a harmonic surface current is prescribed
//
// The following parameters can be changed:

const int P_INIT = 0;             // Initial polynomial degree. NOTE: The meaning is different from
                                  // standard continuous elements in the space H1. Here, P_INIT refers
                                  // to the maximum poly order of the tangential component, and polynomials
                                  // of degree P_INIT + 1 are present in element interiors. P_INIT = 0
                                  // is for Whitney elements.
const bool ALIGN_MESH = false;    // if ALIGN_MESH == true, curvilinear elements aligned with the
                                  // circular load are used, otherwise one uses a non-aligned mesh.
const double THRESHOLD = 0.3;     // This is a quantitative parameter of the adapt(...) function and
                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;           // Adaptive strategy:
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
const double ERR_STOP = 1e-4;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.


// problem constants
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

// boundary conditions
int e_bc_types(int marker)
{
  if (marker == 2) return BC_ESSENTIAL; // perfect conductor BC
  else return BC_NATURAL;               // impedance BC
}

// defining load geometry
bool in_load(double x, double y)
{
  double cx = -0.152994121;
  double cy =  0.030598824;
  double r = 0.043273273;
  if (sqr(cx - x) + sqr(cy - y) < sqr(r)) return true;
  else return false;
}

// gamma as a function of x, y
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


// relative permittivity as a function of x, y
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


// bilinear form
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  cplx ikappa = cplx(0.0, kappa);
  return 1.0/mu_r * int_curl_e_curl_f<Real, Scalar>(n, wt, u, v) -
         ikappa * sqrt(mu_0 / e_0) * int_F_e_f<Real, Scalar>(n, wt, gam, u, v, e) -
         sqr(kappa) * int_F_e_f<Real, Scalar>(n, wt, er, u, v, e);
}

// surface linear form
template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  cplx ii = cplx(0.0, 1.0);
  return ii * omega * J * int_v1<Real, Scalar>(n, wt, v); // just second component of v, since J = (0, J)
}

// error calculation
template<typename Real, typename Scalar>
Scalar hcurl_form_kappa(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_curl_e_curl_f<Scalar, Scalar>(n, wt, u, v) + sqr(kappa) * int_e_f<Scalar, Scalar>(n, wt, u, v);
}

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  if (ALIGN_MESH) mloader.load("oven_load_circle.mesh", &mesh);
  else mloader.load("oven_load_square.mesh", &mesh);

  //mesh.refine_all_elements();

  // initialize the shapeset and the cache
  HcurlShapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create finite element space
  HcurlSpace space(&mesh, &shapeset);
  space.set_bc_types(e_bc_types);
  space.set_uniform_order(P_INIT);

  // enumerate degrees of freedom
  int ndof = assign_dofs(&space);

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form));
  wf.add_liform_surf(0, callback(linear_form_surf));

  // visualize solution and mesh
  VectorView eview("Electric field",0,0,800, 590);
  OrderView ord("Order", 800, 0, 700, 590);

  /*
  // view the basis functions
  VectorBaseView bview;
  vbview.show(&space);
  vbview.wait_for_keypress();
  */

  // matrix solver
  UmfpackSolver solver;

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof_est, graph_cpu_est;

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

    // coarse problem
    LinSystem sys(&wf, &solver);
    sys.set_spaces(1, &space);
    sys.set_pss(1, &pss);
    sys.assemble();
    sys.solve(1, &sln_coarse);

    // time measurement
    cpu_time.tick();

    // show real part of the solution
    AbsFilter abs(&sln_coarse);
    eview.set_min_max_range(0, 4e3);
    eview.show(&abs);
    ord.show(&space);

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // solve the fine mesh problem
    RefSystem ref(&sys);
    ref.assemble();
    ref.solve(1, &sln_fine);

    // calculate error estimate wrt. fine mesh solution
    HcurlOrthoHP hp(1, &space);
    hp.set_biform(0, 0, callback(hcurl_form_kappa));
    double err_est_adapt = hp.calc_error(&sln_coarse, &sln_fine) * 100;
    double err_est_hcurl = hcurl_error(&sln_coarse, &sln_fine) * 100;

    // time measurement
    cpu_time.tick();

    // report results
    info("Error estimate (adapt): %g%%", err_est_adapt);
    info("Error estimate (hcurl): %g%%", err_est_hcurl);

    // add entries to DOF convergence graphs
    graph_dof_est.add_values(space.get_num_dofs(), err_est_hcurl);
    graph_dof_est.save("conv_dof_est.dat");

    // add entries to CPU convergence graphs
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
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // wait for keyboard or mouse input
  View::wait();
  return 0;
}

