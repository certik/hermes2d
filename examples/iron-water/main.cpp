#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

using namespace RefinementSelectors;

//  This example is a standard nuclear engineering benchmark describing an external-force-driven
//  configuration without fissile materials present, using one-group neutron diffusion approximation.
//  It is very similar to example "saphir", the main difference being that the mesh is loaded in
//  the ExodusII format (created for example by Cubit).
//
//  PDE: -div(D(x,y)grad\Phi) + \Sigma_a(x,y)\Phi = Q_{ext}(x,y)
//  where D(x, y) is the diffusion coefficient, \Sigma_a(x,y) the absorption cross-section,
//  and Q_{ext}(x,y) external sources
//
//  Domain: square (0, L)x(0, L) where L = 30c (see mesh file domain.mesh)
//
//  BC:  Zero Dirichlet for the right and top edges ("vacuum boundary")
//       Zero Neumann for the left and bottom edges ("reflection boundary")
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 0;       // Number of initial uniform mesh refinements
const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.6;     // This is a quantitative parameter of the adapt(...) function and
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
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;   // Maximum allowed level of hanging nodes:
                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                  // Note that regular meshes are not supported, this is due to
                                  // their notoriously bad performance.
const double ERR_STOP = 1e-3;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const double CONV_EXP = 1.0;      // Default value is 1.0. This parameter influences the selection of
                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int NDOF_STOP = 60000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// Problem parameters.
double L = 30;                        // edge of square
double L0 = 0.75*0.5*L;               // end of first water layer
double L1 = 0.5*L;                    // end of second water layer
double L2 = 0.75*L;                   // end of iron layer
double Q_EXT = 1.0;                   // neutron source (nonzero in domain 1 only)
double SIGMA_T_WATER = 3.33;          // total cross-section
double SIGMA_T_IRON = 1.33;
double C_WATER = 0.994;               // scattering ratio
double C_IRON = 0.831;
double D_WATER = 1./(3.*SIGMA_T_WATER);  // diffusion coefficient
double D_IRON = 1./(3.*SIGMA_T_IRON);
double SIGMA_A_WATER = SIGMA_T_WATER - C_WATER*SIGMA_T_WATER;  // absorbing cross-section
double SIGMA_A_IRON = SIGMA_T_IRON - C_IRON*SIGMA_T_IRON;

// Materials:
const int WATER_1 = 1;
const int WATER_2 = 2;
const int IRON = 3;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == 1) return BC_NATURAL;
  else return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0.0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  ExodusIIReader mloader;
  if (!mloader.load("iron-water.e", &mesh)) error("ExodusII mesh load failed.");

  // Perform initial uniform mesh refinement.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize the shapeset and the cache.
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // Create finite element space.
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_essential_bc_values(essential_bc_values);
  space.set_uniform_order(P_INIT);

  // Enumerate basis functions.
  int ndof = assign_dofs(&space);

  // Initialize the weak formulation
  WeakForm wf;
  wf.add_biform(bilinear_form_water, bilinear_form_ord, H2D_SYM, WATER_1);
  wf.add_biform(bilinear_form_water, bilinear_form_ord, H2D_SYM, WATER_2);
  wf.add_biform(bilinear_form_iron, bilinear_form_ord, H2D_SYM, IRON);
  wf.add_liform(linear_form_source, linear_form_ord, WATER_1);

  // Visualize solution and mesh.
  ScalarView sview("Coarse solution", 0, 0, 400, 350);
  OrderView  oview("Polynomial orders", 410, 0, 400, 350);

  // Matrix solver.
  UmfpackSolver solver;

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER, &shapeset);

  // Adaptivity loop:
  int as = 1; bool done = false;
  Solution sln_coarse, sln_fine;
  do
    {
    info("---- Adaptivity step %d:", as);

    // Solve the coarse mesh problem.
    LinSystem ls(&wf, &solver);
    ls.set_space(&space);
    ls.set_pss(&pss);
    ls.assemble();
    ls.solve(&sln_coarse);

    // Time measurement.
    cpu_time.tick();

    // View the solution and mesh.
    sview.show(&sln_coarse);
    oview.show(&space);

    // Skip visualization time. 
    cpu_time.tick(H2D_SKIP);

    // Solve the fine mesh problem.
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(&sln_fine);

    // Calculate error estimate wrt. fine mesh solution.
    H1Adapt hp(&space);
    hp.set_solutions(&sln_coarse, &sln_fine);
    double err_est = hp.calc_error() * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est: %g%%", 
         space.get_num_dofs(), rs.get_space(0)->get_num_dofs(), err_est);

    // Add entry to DOF convergence graph.
    graph_dof.add_values(space.get_num_dofs(), err_est);
    graph_dof.save("conv_dof.dat");

    // Add entry to CPU convergence graph.
    graph_cpu.add_values(cpu_time.accumulated(), err_est);
    graph_cpu.save("conv_cpu.dat");

    // If err_est too large, adapt the mesh.
    if (err_est < ERR_STOP) done = true;
    else {
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      ndof = assign_dofs(&space);
      if (ndof >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the fine solution - the final result.
  sview.set_title("Final solution");
  sview.show(&sln_fine);
  oview.set_title("Final orders");
  oview.show(&space);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
