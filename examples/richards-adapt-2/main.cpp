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
const double TAU = 1e-5;                   // Time step.
const double T_FINAL = 5.0;                // Time interval length.

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
const double NEWTON_TOL_COARSE = 0.01;     // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.05;       // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 150;            // Maximum allowed number of Newton iterations.

// Problem parameters.
//double K_S_1 = 0.789; 
//double K_S_2 = 0.469; 
//double K_S_3 = 1e-2; 
//double K_S_4 = 41.143; 
//double K_S_4 = 0.8143; 
double K_S_1 = 0.108;
double K_S_3 = 0.0048;
double K_S_2 = 0.0168;
//double K_S_4 = 41.143; 
double K_S_4 = 1.061;


// double K_S_1 = 1.0; 
// double K_S_2 = 1.0; 
// double K_S_3 = 1.0; 
// //double K_S_4 = 41.143; 
// double K_S_4 = 1.0; 

//double ALPHA_1 = 0.05;
//double ALPHA_2 = 0.05;
//double ALPHA_3 = 0.05;
//double ALPHA_4 = 0.05;


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

double Q_CONST = 0.02;
double STORATIVITY = 0.05;

// Global variables for forms.
double K_S, ALPHA, THETA_R, THETA_S, N, M;

// Material properties.
bool is_in_mat_1(double x, double y) {
  if (y >= 6.0) return true;
  else return false; 
}

bool is_in_mat_2(double x, double y) {
  if (y >= 5.5 && y < 6.0) return true;
  else return false; 
}

bool is_in_mat_4(double x, double y) {
  if (x >= 1.0 && x <= 3.0 && y >= 4.0 && y < 5.0) return true;
  else return false; 
}

bool is_in_mat_3(double x, double y) {
  if (!is_in_mat_1(x, y) && !is_in_mat_2(x, y) && !is_in_mat_4(x, y)) return true;
  else return false; 
}

// Current time.
const double TIME_INIT = 0;                // Initial time.
double TIME = TIME_INIT;

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
  dy = 0;
  return -3.0;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 2.2 - y;
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
  //selector.set_error_weights(2.0, 1.0, sqrt(2.0));
  //selector.set_error_weights(1.0, 1.0, 1.0);

  // Adapt mesh to represent initial condition with given accuracy.
  int proj_norm = 1;  // H1 norm.
  bool verbose = true; 
  bool use_projection = true;
  bool visualization = true; 
  adapt_to_exact_function(&space, init_cond, &selector, THRESHOLD, STRATEGY, 
                          MESH_REGULARITY, ERR_STOP, NDOF_STOP, proj_norm, 
	                  use_projection, verbose, visualization);



  // Solutions for the time stepping and the Newton's method.
  Solution u_prev_time, u_prev_newton;

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(jac_form_vol, jac_form_vol_ord, H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&u_prev_newton, &u_prev_time));
  wf.add_matrix_form_surf(jac_form_surf_1, jac_form_surf_1_ord, BDY_1, Tuple<MeshFunction*>(&u_prev_newton));
  wf.add_matrix_form_surf(jac_form_surf_4, jac_form_surf_4_ord, BDY_4, Tuple<MeshFunction*>(&u_prev_newton));
  wf.add_matrix_form_surf(jac_form_surf_6, jac_form_surf_6_ord, BDY_6, Tuple<MeshFunction*>(&u_prev_newton));
  wf.add_vector_form(res_form_vol, res_form_vol_ord, H2D_ANY, Tuple<MeshFunction*>(&u_prev_newton, &u_prev_time));
  wf.add_vector_form_surf(res_form_surf_1, res_form_surf_1_ord, BDY_1, Tuple<MeshFunction*>(&u_prev_newton));
  wf.add_vector_form_surf(res_form_surf_4, res_form_surf_4_ord, BDY_4, Tuple<MeshFunction*>(&u_prev_newton));
  wf.add_vector_form_surf(res_form_surf_6, res_form_surf_6_ord, BDY_6, Tuple<MeshFunction*>(&u_prev_newton));
  
  // Initialize the nonlinear system.
  NonlinSystem nls(&wf, &space);

  // Error estimate and discrete problem size as a function of physical time.
  SimpleGraph graph_time_err_est, graph_time_dof_est;

  // Set the Dirichlet lift to be the initial solution.
  // The initial vector for the Newton's method will be zero. 
  info("Setting initial vector for the Newton's method zero.");
  u_prev_time.set_dirichlet_lift(&space);             // u_prev_time set equal to init_cond().

  nls.project_global(&u_prev_time, &u_prev_newton);   // Initial vector calculated here.

  // View the projection of the initial condition.
//   ScalarView view("Projection of initial condition", 0, 0, 410, 300);
//   OrderView ordview("Initial mesh", 420, 0, 350, 300);
//   view.fix_scale_width(80);
//   view.show(&u_prev_newton);
//   ordview.show(&space);
  //View::wait(H2DV_WAIT_KEYPRESS);

  // Newton's loop on the coarse mesh.
  info("Solving on coarse mesh.");
  bool verbose2 = true; // Default is false.
  if (!nls.solve_newton(&u_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose2))
    error("Newton's method did not converge.");

  // Store the result in sln_coarse.
  Solution sln_coarse, sln_fine;
  sln_coarse.copy(&u_prev_newton);

  // Time stepping loop.
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    // Updating current time.
    TIME = ts*TAU;

    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);

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

      // Newton's method on fine mesh
      info("Solving on fine mesh.");
      if (!rnls.solve_newton(&u_prev_newton, NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose2))
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

    // write solution into file
    // create filename
    char* filenamecoarse = new char[100];
//    char* filenamefine = new char[100];
//    char* filefinemesh = new char[100];
 //   char* fielcoarsemesh = new char[100];
    sprintf(filenamecoarse, "coarse_%g.dat", TIME);
//    sprintf(filenamefine, "fine_%g.dat", TIME);


 //   sprintf(filefinemesh, "mesh_%g.dat", TIME);

    bool compress = false ;
    sln_coarse.save(  filenamecoarse,  compress );
 //   sln_fine.save(  filenamefine , compress );

  //  mesh.save( filefinemesh ) ;

    
    
    // Visualize the solution and mesh.
    char title[100];
    sprintf(title, "Solution, time level %d", ts);
//     view.set_title(title);
//     view.show(&sln_coarse);
    sprintf(title, "Mesh, time level %d", ts);
//     ordview.set_title(title);
//     ordview.show(&space);

    // Add entries to convergence graphs.
    graph_time_err_est.add_values(ts*TAU, space_err_est);
    graph_time_err_est.save("time_error_est.dat");
    graph_time_dof_est.add_values(ts*TAU, nls.get_num_dofs());
    graph_time_dof_est.save("time_dof_est.dat");

    // Copy new time level solution into u_prev_time.
    u_prev_time.copy(&sln_fine);
  }

  // Wait for all views to be closed.
//   View::wait();
  return 0;
}
