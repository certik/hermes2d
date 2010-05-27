#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

using namespace RefinementSelectors;

// This example employs the multimesh adaptive hp-FEM for linear
// elasticity equations. The domain contains two horizontal
// cracks causing strong singularities at their corners. Each
// displacement component is approximated on an individual mesh.
//
// PDE: Lame equations of linear elasticity.
//
// BC: u_1 = u_2 = 0 on Gamma_1 (left edge)
//     du_2/dn = f on Gamma_2 (upper edge)
//     du_1/dn = du_2/dn = 0 elsewhere, including two horizontal
//               cracks inside the domain. The width of the cracks
//               is currently zero, it can be set in the mesh file
//               via the parameter 'w'.
//
// The following parameters can be changed:

const int INIT_REF_NUM = 0;          // Number of initial uniform mesh refinements.
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
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;      // Maximum allowed level of hanging nodes:
                                     // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                     // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                     // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                     // Note that regular meshes are not supported, this is due to
                                     // their notoriously bad performance.
const double CONV_EXP = 1.0;         // Default value is 1.0. This parameter influences the selection of
                                     // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.5;         // Stopping criterion for adaptivity (rel. error tolerance between the
                                     // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;         // Adaptivity process stops when the number of degrees of freedom grows.

// Problem parameters.
const double E  = 200e9;             // Young modulus for steel: 200 GPa.
const double nu = 0.3;               // Poisson ratio.
const double f  = 1e3;               // Load force.
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));

// Boundary markers.
const int BDY_LEFT = 1;
const int BDY_TOP = 2;

// Boundary condition types.
BCType bc_types_xy(int marker)
  { return (marker == BDY_LEFT) ? BC_ESSENTIAL : BC_NATURAL; }

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh xmesh, ymesh;
  H2DReader mloader;
  mloader.load("crack.mesh", &xmesh);

  // Perform initial uniform mesh refinement.
  for (int i=0; i < INIT_REF_NUM; i++) xmesh.refine_all_elements();

  // Create initial mesh for the vertical displacement component.
  // This also initializes the multimesh hp-FEM.
  ymesh.copy(&xmesh);

  // Initialize the shapeset and the cache.
  H1Shapeset shapeset;
  PrecalcShapeset xpss(&shapeset);
  PrecalcShapeset ypss(&shapeset);

  // Create the x displacement space.
  H1Space xdisp(&xmesh, &shapeset);
  xdisp.set_bc_types(bc_types_xy);
  xdisp.set_uniform_order(P_INIT);

  // Create the y displacement space.
  H1Space ydisp(MULTI ? &ymesh : &xmesh, &shapeset);
  ydisp.set_bc_types(bc_types_xy);
  ydisp.set_uniform_order(P_INIT);

  // Enumerate degrees of freedom.
  int ndof = assign_dofs(2, &xdisp, &ydisp);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_biform(0, 0, callback(bilinear_form_0_0), H2D_SYM);
  wf.add_biform(0, 1, callback(bilinear_form_0_1), H2D_SYM);
  wf.add_biform(1, 1, callback(bilinear_form_1_1), H2D_SYM);
  wf.add_liform_surf(1, callback(linear_form_surf_1), BDY_TOP);

  // Initialize views.
  ScalarView sview("Von Mises stress [Pa]", 0, 355, 900, 300);
  OrderView  xoview("X polynomial orders", 0, 0, 900, 300);
  OrderView  yoview("Y polynomial orders", 910, 0, 900, 300);

  // Matrix solver.
  UmfpackSolver solver;

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER, &shapeset);

  // Adaptivity loop:
  int as = 1; bool done = false;
  Solution sln_x_coarse, sln_y_coarse, sln_x_fine, sln_y_fine;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Solve the coarse mesh problem.
    LinSystem ls(&wf, &solver);
    ls.set_spaces(2, &xdisp, &ydisp);
    ls.set_pss(2, &xpss, &ypss);
    ls.assemble();
    ls.solve(2, &sln_x_coarse, &sln_y_coarse);

    // Time measurement.
    cpu_time.tick();

    // Visualize the solution and meshes.
    VonMisesFilter stress(&sln_x_coarse, &sln_y_coarse, mu, lambda);
    //sview.set_min_max_range(0, 3e4);
    sview.show(&stress, H2D_EPS_HIGH);
    xoview.show(&xdisp);
    yoview.show(&ydisp);

    // Time measurement.
    cpu_time.tick(H2D_SKIP);

    // Solve the fine mesh problem.
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(2, &sln_x_fine, &sln_y_fine);

    // Calculate error estimate wrt. fine mesh solution in energy norm.
    H1Adapt hp(Tuple<Space*>(&xdisp, &ydisp));
    hp.set_solutions(Tuple<Solution*>(&sln_x_coarse, &sln_y_coarse), Tuple<Solution*>(&sln_x_fine, &sln_y_fine));
    hp.set_biform(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
    hp.set_biform(0, 1, bilinear_form_0_1<scalar, scalar>, bilinear_form_0_1<Ord, Ord>);
    hp.set_biform(1, 0, bilinear_form_1_0<scalar, scalar>, bilinear_form_1_0<Ord, Ord>);
    hp.set_biform(1, 1, bilinear_form_1_1<scalar, scalar>, bilinear_form_1_1<Ord, Ord>);
    double err_est = hp.calc_error(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;

    // Time measurement.
    cpu_time.tick();

    // Report results.
    info("ndof_x_coarse: %d, ndof_x_fine: %d", 
         xdisp.get_num_dofs(), rs.get_space(0)->get_num_dofs());
    info("ndof_y_coarse: %d, ndof_y_fine: %d", 
         ydisp.get_num_dofs(), rs.get_space(1)->get_num_dofs());
    info("err_est: %g%%", err_est);

    // Add entry to DOF convergence graph.
    graph_dof.add_values(xdisp.get_num_dofs() + ydisp.get_num_dofs(), err_est);
    graph_dof.save("conv_dof.dat");

    // Add entry to CPU convergence graph.
    graph_cpu.add_values(cpu_time.accumulated(), err_est);
    graph_cpu.save("conv_cpu.dat");

    // If err_est too large, adapt the mesh.
    if (err_est < ERR_STOP || xdisp.get_num_dofs() + ydisp.get_num_dofs() >= NDOF_STOP) done = true;
    else {
      done = hp.adapt(&selector, MULTI ? THRESHOLD_MULTI : THRESHOLD_SINGLE, STRATEGY, MESH_REGULARITY, SAME_ORDERS);
      ndof = assign_dofs(2, &xdisp, &ydisp);
      if (ndof >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (!done);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
