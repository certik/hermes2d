#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

// This example employs the multimesh adaptive hp-FEM for linear
// elasticity equations. The domain contains two horizontal
// cracks causing strong singularities at their corners. Each
// displacement component is approximated on an individual mesh.
//
// PDE: Lame equations of linear elasticity
//
// BC: u_1 = u_2 = 0 on Gamma_1 (left edge)
//     du_2/dn = f on Gamma_2 (upper edge)
//     du_1/dn = du_2/dn = 0 elsewhere, including two horizontal
//               cracks inside the domain. The width of the cracks
//               is currently zero, it can be set in the mesh file
//               via the parameter 'w'.
//
// The following parameters can be changed:

const int P_INIT = 2;                // Initial polynomial degree of all mesh elements.
const bool MULTI = true;             // true = use multi-mesh, false = use single-mesh.
                                     // Note: in the single mesh option, the meshes are
                                     // forced to be geometrically the same but the
                                     // polynomial degrees can still vary.
const bool SAME_ORDERS = false;      // true = when single mesh is used it forces same pol.
                                     // orders for components
                                     // when multi mesh used, parameter is ignored
const double THRESHOLD_MULTI = 0.35; // error threshold for element refinement (multi-mesh)
const double THRESHOLD_SINGLE = 0.7; // error threshold for element refinement (single-mesh)
const int STRATEGY = 0;              // Adaptive strategy:
                                     // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                     //   error is processed. If more elements have similar errors, refine
                                     //   all to keep the mesh symmetric.
                                     // STRATEGY = 1 ... refine all elements whose error is larger
                                     //   than THRESHOLD times maximum element error.
                                     // STRATEGY = 2 ... refine all elements whose error is larger
                                     //   than THRESHOLD.
                                     // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int ADAPT_TYPE = 0;            // Type of automatic adaptivity:
                                     // ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                     // ADAPT_TYPE = 1 ... adaptive h-FEM,
                                     // ADAPT_TYPE = 2 ... adaptive p-FEM.
const bool ISO_ONLY = false;         // Isotropic refinement flag (concerns quadrilateral elements only).
                                     // ISO_ONLY = false ... anisotropic refinement of quad elements
                                     // is allowed (default),
                                     // ISO_ONLY = true ... only isotropic refinements of quad elements
                                     // are allowed.
const int MESH_REGULARITY = -1;      // Maximum allowed level of hanging nodes:
                                     // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                     // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                     // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                     // Note that regular meshes are not supported, this is due to
                                     // their notoriously bad performance.
const double CONV_EXP = 1.0;         // Default value is 1.0. This parameter influences the selection of
                                     // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int MAX_ORDER = 10;            // Maximum polynomial order used during adaptivity.
const double ERR_STOP = 1e-2;        // Stopping criterion for adaptivity (rel. error tolerance between the
                                     // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;         // Adaptivity process stops when the number of degrees of freedom grows

// problem constants
const double E  = 200e9;  // Young modulus for steel: 200 GPa
const double nu = 0.3;    // Poisson ratio
const double f  = 1e3;    // load force
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));

// boundary markers
const int marker_left = 1;
const int marker_top = 2;

// boundary conditions
int bc_types_xy(int marker)
  { return (marker == marker_left) ? BC_ESSENTIAL : BC_NATURAL; }

// linear and bilinear forms
template<typename Real, typename Scalar>
Scalar bilinear_form_0_0(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (lambda + 2*mu) * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
                      mu * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_0_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return lambda * int_dudy_dvdx<Real, Scalar>(n, wt, u, v) +
             mu * int_dudx_dvdy<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_0(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return     mu * int_dudy_dvdx<Real, Scalar>(n, wt, u, v) +
         lambda * int_dudx_dvdy<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return              mu * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
         (lambda + 2*mu) * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_1(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return -f * int_v<Real, Scalar>(n, wt, v);
}

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh xmesh, ymesh;
  H2DReader mloader;
  mloader.load("crack.mesh", &xmesh);
  ymesh.copy(&xmesh);          // this defines the common master mesh for
                               // both displacement fields

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset xpss(&shapeset);
  PrecalcShapeset ypss(&shapeset);

  // create the x displacement space
  H1Space xdisp(&xmesh, &shapeset);
  xdisp.set_bc_types(bc_types_xy);
  xdisp.set_uniform_order(P_INIT);

  // create the y displacement space
  H1Space ydisp(MULTI ? &ymesh : &xmesh, &shapeset);
  ydisp.set_bc_types(bc_types_xy);
  ydisp.set_uniform_order(P_INIT);

  // enumerate basis functions
  int ndof = assign_dofs(2, &xdisp, &ydisp);
  printf("xdof=%d, ydof=%d\n", xdisp.get_num_dofs(), ydisp.get_num_dofs());

  // initialize the weak formulation
  WeakForm wf(2);
  wf.add_biform(0, 0, callback(bilinear_form_0_0), H2D_SYM);
  wf.add_biform(0, 1, callback(bilinear_form_0_1), H2D_SYM);
  wf.add_biform(1, 1, callback(bilinear_form_1_1), H2D_SYM);
  wf.add_liform_surf(1, callback(linear_form_surf_1), marker_top);

  // visualize solution and mesh
  ScalarView sview("Von Mises stress [Pa]", 0, 355, 900, 300);
  OrderView  xoview("X polynomial orders", 0, 0, 900, 300);
  OrderView  yoview("Y polynomial orders", 910, 0, 900, 300);

  // matrix solver
  UmfpackSolver solver;

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof, graph_cpu;

  // adaptivity loop
  int it = 1;
  bool done = false;
  TimePeriod cpu_time;
  Solution sln_x_coarse, sln_y_coarse, sln_x_fine, sln_y_fine;
  do
  {
    info("!---- Adaptivity step %d ---------------------------------------------", it); it++;

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // solve the coarse mesh problem
    LinSystem ls(&wf, &solver);
    ls.set_spaces(2, &xdisp, &ydisp);
    ls.set_pss(2, &xpss, &ypss);
    ls.assemble();
    ls.solve(2, &sln_x_coarse, &sln_y_coarse);

    // time measurement
    cpu_time.tick();

    // visualize the solution
    VonMisesFilter stress(&sln_x_coarse, &sln_y_coarse, mu, lambda);
    //sview.set_min_max_range(0, 3e4);
    sview.show(&stress, H2D_EPS_HIGH);
    xoview.show(&xdisp);
    yoview.show(&ydisp);

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // solve the fine (reference) problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(2, &sln_x_fine, &sln_y_fine);

    // calculate error estimate wrt. fine mesh solution in energy norm
    H1OrthoHP hp(2, &xdisp, &ydisp);
    hp.set_biform(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
    hp.set_biform(0, 1, bilinear_form_0_1<scalar, scalar>, bilinear_form_0_1<Ord, Ord>);
    hp.set_biform(1, 0, bilinear_form_1_0<scalar, scalar>, bilinear_form_1_0<Ord, Ord>);
    hp.set_biform(1, 1, bilinear_form_1_1<scalar, scalar>, bilinear_form_1_1<Ord, Ord>);
    double err_est = hp.calc_error_2(&sln_x_coarse, &sln_y_coarse, &sln_x_fine, &sln_y_fine) * 100;

    // time measurement
    cpu_time.tick();

    // report results
    info("Error estimate: %g %%", err_est);

    // add entry to DOF convergence graph
    graph_dof.add_values(xdisp.get_num_dofs() + ydisp.get_num_dofs(), err_est);
    graph_dof.save("conv_dof.dat");

    // add entry to CPU convergence graph
    graph_cpu.add_values(cpu_time.accumulated(), err_est);
    graph_cpu.save("conv_cpu.dat");

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP || xdisp.get_num_dofs() + ydisp.get_num_dofs() >= NDOF_STOP) done = true;
    else {
      hp.adapt(MULTI ? THRESHOLD_MULTI : THRESHOLD_SINGLE, STRATEGY, ADAPT_TYPE, ISO_ONLY, MESH_REGULARITY, CONV_EXP, MAX_ORDER, SAME_ORDERS);
      ndof = assign_dofs(2, &xdisp, &ydisp);
      printf("xdof=%d, ydof=%d\n", xdisp.get_num_dofs(), ydisp.get_num_dofs());
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
