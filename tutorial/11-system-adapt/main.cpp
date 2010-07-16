#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This example explains how to use the multimesh adaptive hp-FEM,
// where different physical fields (or solution components) can be
// approximated using different meshes and equipped with mutually
// independent adaptivity mechanisms. For the tutorial purposes,
// we manufactured an exact solution for a simplified version of
// the FitzHugh-Nagumo equation. This equation, in its full form,
// is a prominent example of activator-inhibitor systems in two-component
// reaction-diffusion equations, It describes a prototype of an
// excitable system (e.g., a neuron).
//
// PDE: Linearized FitzHugh-Nagumo equation
//      -d_u^2 \Delta u - f(u) + \sigma v = g_1,
//      -d_v^2 \Delta v - u + v = g_2.
// In the original equation, f(u) = \lambda u - u^3 - \kappa. For
// simplicity, here we just take f(u) = u.
//
// Domain: Square (-1,1)^2.
//
// BC: Both solution components are zero on the boundary.
//
// Exact solution: The functions g_1 and g_2 were calculated so that
//                 the exact solution is:
//        u(x,y) = U(x)*U(y) where U(t) = cos(M_PI*t/2)
//        v(x,y) = V(x)V(y) where V(t) = 1 - (exp(K*t)+exp(-K*t))/(exp(K) + exp(-K))
// Note: V(t) is the exact solution of the 1D singularly perturbed equation
//       -u'' + K*K*u = K*K in (-1, 1) with zero Dirichlet BC.
//
// The following parameters can be changed: In particular, compare hp- and
// h-adaptivity via the CAND_LIST option, and compare the multi-mesh vs.
// single-mesh using the MULTI parameter.

const int P_INIT_U = 2;                  // Initial polynomial degree for u.
const int P_INIT_V = 2;                  // Initial polynomial degree for v.
const int INIT_REF_BDY = 3;              // Number of initial boundary refinements
const bool MULTI = true;                 // MULTI = true  ... use multi-mesh,
                                         // MULTI = false ... use single-mesh.
                                         // Note: In the single mesh option, the meshes are
                                         // forced to be geometrically the same but the
                                         // polynomial degrees can still vary.
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
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
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See the User Documentation for details.
const int MESH_REGULARITY = -1;  // Maximum allowed level of hanging nodes:
                                 // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                 // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                 // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                 // Note that regular meshes are not supported, this is due to
                                 // their notoriously bad performance.
const double CONV_EXP = 1;       // Default value is 1.0. This parameter influences the selection of
                                 // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int MAX_ORDER = 10;        // Maximum allowed element degree
const double ERR_STOP = 0.5;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                 // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;     // Adaptivity process stops when the number of degrees of freedom grows over
                                 // this limit. This is mainly to prevent h-adaptivity to go on forever.

// Problem parameters.
const double D_u = 1;
const double D_v = 1;
const double SIGMA = 1;
const double LAMBDA = 1;
const double KAPPA = 1;
const double K = 100;

// Boundary condition types.
BCType bc_types(int marker) { return BC_ESSENTIAL; }

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y) { return 0;}

// Exact solution.
#include "exact_solution.cpp"

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh u_mesh, v_mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &u_mesh);
  if (MULTI == false) u_mesh.refine_towards_boundary(1, INIT_REF_BDY);

  // Create initial mesh (master mesh).
  v_mesh.copy(&u_mesh);

  // Initial mesh refinements in the v_mesh towards the boundary.
  if (MULTI == true) v_mesh.refine_towards_boundary(1, INIT_REF_BDY);

  // Create H1 spaces with default shapeset for both displacement components.
  H1Space u_space(&u_mesh, bc_types, essential_bc_values, P_INIT_U);
  H1Space v_space(MULTI ? &v_mesh : &u_mesh, bc_types, essential_bc_values, P_INIT_V);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0));
  wf.add_matrix_form(0, 1, callback(bilinear_form_0_1));
  wf.add_matrix_form(1, 0, callback(bilinear_form_1_0));
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1));
  wf.add_vector_form(0, linear_form_0, linear_form_0_ord);
  wf.add_vector_form(1, linear_form_1, linear_form_1_ord);

  // Initialize views.
  OrderView  uoview("Coarse mesh for u", 0, 0, 360, 300);
  OrderView  voview("Coarse mesh for v", 370, 0, 360, 300);
  ScalarView uview("Coarse mesh solution u", 740, 0, 400, 300);
  ScalarView vview("Coarse mesh solution v", 1150, 0, 400, 300);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;  
  init_matrix_solver(SOLVER_UMFPACK, u_space.get_num_dofs() + v_space.get_num_dofs(), 
                     mat, rhs, solver);

  // Adaptivity loop.
  Solution u_sln, v_sln;
  Solution ref_u_sln, ref_v_sln;
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);
    info("Solving on reference meshes.");

    // Construct globally refined reference meshes.
    Mesh ref_u_mesh, ref_v_mesh;
    ref_u_mesh.copy(&u_mesh);
    ref_v_mesh.copy(&v_mesh);
    ref_u_mesh.refine_all_elements();
    ref_v_mesh.refine_all_elements();

    // Setup spaces for the reference solution.
    Space *ref_u_space = u_space.dup(&ref_u_mesh);
    Space *ref_v_space = v_space.dup(&ref_v_mesh);
    int order_increase = 1;
    ref_u_space->copy_orders(&u_space, order_increase);
    ref_v_space->copy_orders(&v_space, order_increase);
 
    // Solve the reference problem.
    solve_linear(Tuple<Space *>(ref_u_space, ref_v_space), &wf, 
                 Tuple<Solution *>(&ref_u_sln, &ref_v_sln), SOLVER_UMFPACK);

    // Project the reference solutions on the coarse meshes.
    info("Projecting reference solutions on coarse meshes.");
    project_global(Tuple<Space *>(&u_space, &v_space), 
                   Tuple<MeshFunction*>(&ref_u_sln, &ref_v_sln), 
                   Tuple<Solution*>(&u_sln, &v_sln));

    // Time measurement.
    cpu_time.tick();

    // View the solutions and meshes.
    info("u_dof_coarse: %d, v_dof_coarse: %d", u_space.get_num_dofs(), v_space.get_num_dofs());
    info("ref_u_dof: %d, ref_v_dof: %d", ref_u_space->get_num_dofs(), ref_v_space->get_num_dofs());
    uview.show(&u_sln);
    vview.show(&v_sln);
    uoview.show(&u_space);
    voview.show(&v_space);

    // Time measurement.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors and total error estimate.
    info("Calculating error (est).");
    H1Adapt hp(Tuple<Space *>(&u_space, &v_space));
    hp.set_solutions(Tuple<Solution*>(&u_sln, &v_sln), 
                     Tuple<Solution*>(&ref_u_sln, &ref_v_sln));
    hp.set_error_form(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
    hp.set_error_form(0, 1, bilinear_form_0_1<scalar, scalar>, bilinear_form_0_1<Ord, Ord>);
    hp.set_error_form(1, 0, bilinear_form_1_0<scalar, scalar>, bilinear_form_1_0<Ord, Ord>);
    hp.set_error_form(1, 1, bilinear_form_1_1<scalar, scalar>, bilinear_form_1_1<Ord, Ord>);
    double err_est = hp.calc_error(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_ABS) * 100;

    // Calculate error wrt. exact solution.
    info("Calculating error (exact).");
    ExactSolution uexact(&u_mesh, u_exact);
    ExactSolution vexact(&v_mesh, v_exact);
    double u_error = h1_error(&u_sln, &uexact) * 100;
    double v_error = h1_error(&v_sln, &vexact) * 100;
    double error = std::max(u_error, v_error);

    // Report results.
    info("Exact solution error for u (H1 norm): %g%%", u_error);
    info("Exact solution error for v (H1 norm): %g%%", v_error);
    info("Exact solution error (maximum): %g%%", error);
    info("Estimate of error wrt. ref. solution (energy norm): %g%%", err_est);

    // Add entry to DOF convergence graph.
    graph_dof.add_values(u_space.get_num_dofs() + v_space.get_num_dofs(), error);
    if (MULTI == true) graph_dof.save("conv_dof_m.dat");
    else graph_dof.save("conv_dof_s.dat");

    // Add entry to CPU convergence graph.
    graph_cpu.add_values(cpu_time.accumulated(), error);
    if (MULTI == true) graph_cpu.save("conv_cpu_m.dat");
    else graph_cpu.save("conv_cpu_s.dat");

    // If err_est too large, adapt the mesh.
    if (error < ERR_STOP) done = true;
    else {
      info("Adapting coarse meshes.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY, MULTI == true ? false : true);
      if (u_space.get_num_dofs() + v_space.get_num_dofs() >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (!done);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the fine solution - the final result.
  uview.set_title("Reference solution u");
  uview.show_mesh(false);
  uview.show(&ref_u_sln);
  vview.set_title("Reference solution v");
  vview.show_mesh(false);
  vview.show(&ref_v_sln);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
