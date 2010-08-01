#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example uses adaptivity with dynamical meshes to solve
//  the time-dependent Richard's equation. The time discretization 
//  is backward Euler, and the Newton's method is applied to solve 
//  the nonlinear problem in each time step. 
//
//  PDE: C(h)dh/dt - div(K(h)grad(h)) - (dK/dh)*(dh/dy) = 0
//  where K(h) = K_S*exp(alpha*h)                          for h < 0,
//        K(h) = K_S                                       for h >= 0,
//        C(h) = alpha*(theta_s - theta_r)*exp(alpha*h)    for h < 0,
//        C(h) = alpha*(theta_s - theta_r)                 for h >= 0.
//
//  Domain: rectangle (0, 8) x (0, 6.5).
//
//  BC: Dirichlet, given by the initial condition.
//  IC: See the function init_cond().
//
//  The following parameters can be changed:

// If this is defined, use van Genuchten's constitutive relations, otherwise use Gardner's.
// #define CONSTITUTIVE_GENUCHTEN

const int P_INIT = 1;                      // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 0;                // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 0;            // Number of initial mesh refinements towards the top edge.

// Adaptivity
const int UNREF_FREQ = 1;                  // Every UNREF_FREQth time step the mesh is unrefined.
const double THRESHOLD = 0.3;              // This is a quantitative parameter of the adapt(...) function and
                                           // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                    // Adaptive strategy:
                                           // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                           //   error is processed. If more elements have similar errors, refine
                                           //   all to keep the mesh symmetric.
                                           // STRATEGY = 1 ... refine all elements whose error is larger
                                           //   than THRESHOLD times maximum element error.
                                           // STRATEGY = 2 ... refine all elements whose error is larger
                                           //   than THRESHOLD.
                                           // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;   // Predefined list of element refinement candidates. Possible values are
                                           // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                           // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                           // See the User Documentation for details.
const int MESH_REGULARITY = -1;            // Maximum allowed level of hanging nodes:
                                           // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                           // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                           // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                           // Note that regular meshes are not supported, this is due to
                                           // their notoriously bad performance.
const double CONV_EXP = 1.0;               // Default value is 1.0. This parameter influences the selection of
                                           // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.5;               // Stopping criterion for adaptivity (rel. error tolerance between the
                                           // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.

// Newton's method
const double NEWTON_TOL_COARSE = 0.0001;   // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.0005;     // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 50;            // Maximum allowed number of Newton iterations.

// Problem parameters.
const double TAU = 5e-3;                   // Time step.
const double STARTUP_TIME = 1.1e-2;        // Start-up time for time-dependent Dirichlet boundary condition.
const double T_FINAL = 5.0;                // Time interval length.
double TIME = 0;                           // Global time variable initialized with first time step.
double H_INIT = -9.5;                      // Initial pressure head.
double H_ELEVATION = 5.2;

double K_S_1 = 0.108;
double K_S_3 = 0.0048;
double K_S_2 = 0.0168;
double K_S_4 = 1.061;

double ALPHA_1 = 0.01;
double ALPHA_3 = 0.005;
double ALPHA_2 = 0.01;
double ALPHA_4 = 0.05;

double THETA_R_1 = 0.1020;
double THETA_R_2 = 0.09849;
double THETA_R_3 = 0.08590;
double THETA_R_4 = 0.08590;

double THETA_S_1 = 0.4570;
double THETA_S_2 = 0.4510;
double THETA_S_3 = 0.4650;
double THETA_S_4 = 0.5650;

double N_1 = 1.982;
double N_2 = 1.632; 
double N_3 = 5.0;
double N_4 = 5.0;

double M_1 = 0.49546;
double M_2 = 0.38726;
double M_3 = 0.8;
double M_4 = 0.8;

double Q_MAX_VALUE = 0.07;         // Maximum value, used in function q_function(); 
double q_function() {
  if (STARTUP_TIME > TIME) return Q_MAX_VALUE * TIME / STARTUP_TIME;
  else return Q_MAX_VALUE;
}

double STORATIVITY = 0.05;

// Global variables for forms.
double K_S, ALPHA, THETA_R, THETA_S, N, M;

// Material properties.
bool is_in_mat_1(double x, double y) {
  if (y >= -0.5) return true;
  else return false; 
}

bool is_in_mat_2(double x, double y) {
  if (y >= -1.0 && y < -0.5) return true;
  else return false; 
}

bool is_in_mat_4(double x, double y) {
  if (x >= 1.0 && x <= 3.0 && y >= -2.5 && y < -1.5) return true;
  else return false; 
}

bool is_in_mat_3(double x, double y) {
  if (!is_in_mat_1(x, y) && !is_in_mat_2(x, y) && !is_in_mat_4(x, y)) return true;
  else return false; 
}

#ifdef CONSTITUTIVE_GENUCHTEN
#include "constitutive_genuchten.cpp"
#else
#include "constitutive_gardner.cpp"
#endif

// Boundary markers.
int BDY_1 = 1;
int BDY_3 = 3;
int BDY_4 = 4;
int BDY_6 = 6;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == BDY_3) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Initial condition.
double init_cond(double x, double y, double& dx, double& dy) {
  dx = 0;
  dy = -1;
  return -y + H_INIT;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  if (STARTUP_TIME > TIME) return -y + H_INIT + TIME/STARTUP_TIME*H_ELEVATION;
  else return -y + H_INIT + H_ELEVATION;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &basemesh);

  // Perform initial mesh refinements.
  mesh.copy(&basemesh);
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(3, INIT_REF_NUM_BDY);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Solutions for the time stepping and the Newton's method.
  Solution u_prev_time, u_prev_newton;

  // Adapt mesh to represent initial condition with given accuracy.
  int proj_norm = 1;  // H1 norm.
  bool verbose0 = true; 
  bool use_projection = true;
  bool visualization = false;
  double err_stop_init_cond = 0.1 * ERR_STOP; 
  adapt_to_exact_function(&space, init_cond, &selector, THRESHOLD, STRATEGY, 
                          MESH_REGULARITY, ERR_STOP, NDOF_STOP, proj_norm, 
	                  use_projection, verbose0, visualization, &u_prev_time);

  // Initialize views.
  ScalarView sview("Solution", 0, 0, 500, 350);
  sview.set_min_max_range(-9.5, 2.2);
  OrderView oview("Mesh", 510, 0, 500, 350);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(jac_form_vol, jac_form_vol_ord, H2D_UNSYM, H2D_ANY, 
                     Tuple<MeshFunction*>(&u_prev_newton, &u_prev_time));
  wf.add_matrix_form_surf(jac_form_surf_1, jac_form_surf_1_ord, BDY_1, 
                     Tuple<MeshFunction*>(&u_prev_newton));
  wf.add_matrix_form_surf(jac_form_surf_4, jac_form_surf_4_ord, BDY_4, 
                     Tuple<MeshFunction*>(&u_prev_newton));
  wf.add_matrix_form_surf(jac_form_surf_6, jac_form_surf_6_ord, BDY_6, 
                     Tuple<MeshFunction*>(&u_prev_newton));
  wf.add_vector_form(res_form_vol, res_form_vol_ord, H2D_ANY, 
                     Tuple<MeshFunction*>(&u_prev_newton, &u_prev_time));
  wf.add_vector_form_surf(res_form_surf_1, res_form_surf_1_ord, BDY_1, 
                     Tuple<MeshFunction*>(&u_prev_newton));
  wf.add_vector_form_surf(res_form_surf_4, res_form_surf_4_ord, BDY_4, 
                     Tuple<MeshFunction*>(&u_prev_newton));
  wf.add_vector_form_surf(res_form_surf_6, res_form_surf_6_ord, BDY_6, 
                     Tuple<MeshFunction*>(&u_prev_newton));
  
  // Initialize the nonlinear system.
  NonlinSystem nls(&wf, &space);

  // Error estimate and discrete problem size as a function of physical time.
  SimpleGraph graph_time_err_est, graph_time_dof_est;

  // Calculating initial vector for Newton.
  info("Projecting initial condition to obtain coefficient vector for Newton on coarse mesh.");
  nls.project_global(&u_prev_time, &u_prev_newton);   // Initial vector calculated here.

  // Newton's loop (one time step) on the coarse mesh.
  info("Solving on coarse mesh.");
  bool verbose = true; // Default is false.
  if (!nls.solve_newton(&u_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
    error("Newton's method did not converge.");

  // Store the result in sln_coarse.
  Solution sln_coarse, sln_fine;
  sln_coarse.copy(&u_prev_newton);

  // Time stepping loop.
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    // Update time-dependent Dirichlet BC values.
    if (TIME <= STARTUP_TIME) space.update_essential_bc_values();

    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
      mesh.refine_towards_boundary(3, INIT_REF_NUM_BDY);
  
      // Project fine mesh solution on the globally derefined mesh.
      info("---- Time step %d:", ts);
      info("Projecting fine mesh solution on globally derefined mesh for error calculation.");
      nls.project_global(&sln_fine, &u_prev_newton);

      // Store the result in sln_coarse.
      sln_coarse.copy(&u_prev_newton);
    }

    // Adaptivity loop (in space):
    bool done = false;
    double space_err_est;
    int as = 1;
    do
    {
      info("---- Time step %d, adaptivity step %d:", ts, as);

      // Initialize reference nonlinear system.
      RefSystem rnls(&nls);

      // Set initial condition for the Newton's method on the fine mesh.
      if (as == 1) {
        info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
        rnls.project_global(&sln_coarse, &u_prev_newton);
      }
      else {
        info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
        rnls.project_global(&sln_fine, &u_prev_newton);
      }

      // Newton's method (one time step) on fine mesh
      info("Solving on fine mesh.");
      if (!rnls.solve_newton(&u_prev_newton, NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose))
        error("Newton's method did not converge.");

      // Store the result in sln_fine.
      sln_fine.copy(&u_prev_newton);

      // Calculate error estimate wrt. fine mesh solution.
      info("Calculating error (est).");
      H1Adapt hp(&nls);
      hp.set_solutions(&sln_coarse, &sln_fine);
      space_err_est = hp.calc_error() * 100;

      info("ndof_coarse: %d, ndof_fine: %d, space_err_est: %g%%", 
	   nls.get_num_dofs(), rnls.get_num_dofs(), space_err_est);

      // If space_err_est too large, adapt the mesh.
      if (space_err_est < ERR_STOP) done = true;
      else {
        info("Adapting coarse mesh.");
        done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (nls.get_num_dofs() >= NDOF_STOP) {
          done = true;
          break;
        }

        // Project the fine mesh solution on the new coarse mesh.
        info("Projecting fine mesh solution on coarse mesh for error calculation.");
        nls.project_global(&sln_fine, &u_prev_newton);

        // Store the result in sln_coarse.
        sln_coarse.copy(&u_prev_newton);

        as++;

      }
    }
    while (!done);

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Solution, t = %g", TIME);
    sview.set_title(title);
    sview.show(&sln_coarse);
    sprintf(title, "Mesh, t = %g", TIME);
    oview.set_title(title);
    oview.show(&space);

    /*
    //Write solution data into file.
    bool compress = false ;
    char* filenamecoarse = new char[100];
    sprintf(filenamecoarse, "coarse_%g.dat", TIME);
    sln_coarse.save(  filenamecoarse,  compress );
    //char* filenamefine = new char[100];
    //sprintf(filenamefine, "fine_%g.dat", TIME);
    //sln_fine.save(filenamefine , compress );
    //char* filefinemesh = new char[100];
    //sprintf(filefinemesh, "mesh_%g.dat", TIME);
    //mesh.save( filefinemesh ) ;
    */

    // Add entries to convergence graphs.
    graph_time_err_est.add_values(ts*TAU, space_err_est);
    graph_time_err_est.save("time_error_est.dat");
    graph_time_dof_est.add_values(ts*TAU, nls.get_num_dofs());
    graph_time_dof_est.save("time_dof_est.dat");

    // Copy new time level solution into u_prev_time.
    u_prev_time.copy(&sln_fine);

    TIME += TAU;
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
