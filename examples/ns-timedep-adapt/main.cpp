#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// The time-dependent laminar incompressible Navier-Stokes equations are
// discretized in time via the implicit Euler method. The Newton's method 
// is used to solve the nonlinear problem at each time step. We show how
// to use discontinuous ($L^2$) elements for pressure and thus make the
// velocity discreetely divergence free. Comparison to approximating the
// pressure with the standard (continuous) Taylor-Hood elements is enabled.
// The Reynolds number Re = 200 which is embarrassingly low. You
// can increase it but then you will need to make the mesh finer, and the
// computation will take more time.
//
// PDE: incompressible Navier-Stokes equations in the form
// \partial v / \partial t - \Delta v / Re + (v \cdot \nabla) v + \nabla p = 0,
// div v = 0
//
// BC: u_1 is a time-dependent constant and u_2 = 0 on Gamma_4 (inlet)
//     u_1 = u_2 = 0 on Gamma_1 (bottom), Gamma_3 (top) and Gamma_5 (obstacle)
//     "do nothing" on Gamma_2 (outlet)
//
// Geometry: Rectangular channel containing an off-axis circular obstacle. The
//           radius and position of the circle, as well as other geometry
//           parameters can be changed in the mesh file "domain.mesh".
//
// The following parameters can be changed:

const bool SOLVE_ON_COARSE_MESH = false; // true... Newton is done on coarse mesh in every adaptivity step.
                                         // false...Newton is done on coarse mesh only once, then projection
                                         // of the fine mesh solution to coarse mesh is used.
const int INIT_REF_NUM = 0;              // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 3;          // Number of initial mesh refinements towards boundary.
#define PRESSURE_IN_L2                   // If this is defined, the pressure is approximated using
                                         // discontinuous L2 elements (making the velocity discreetely
                                         // divergence-free, more accurate than using a continuous
                                         // pressure approximation). Otherwise the standard continuous
                                         // elements are used. The results are striking - check the
                                         // tutorial for comparisons.
const int P_INIT_VEL = 2;                // Initial polynomial degree for velocity components
const int P_INIT_PRESSURE = 1;           // Initial polynomial degree for pressure
                                         // Note: P_INIT_VEL should always be greater than
                                         // P_INIT_PRESSURE because of the inf-sup condition

// Adaptivity
const int UNREF_FREQ = 1;        // Every UNREF_FREQth time step the mesh is unrefined.
const double THRESHOLD = 0.3;    // This is a quantitative parameter of the adapt(...) function and
                                 // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;          // Adaptive strategy:
                                 // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                 //   error is processed. If more elements have similar errors, refine
                                 //   all to keep the mesh symmetric.
                                 // STRATEGY = 1 ... refine all elements whose error is larger
                                 //   than THRESHOLD times maximum element error.
                                 // STRATEGY = 2 ... refine all elements whose error is larger
                                 //   than THRESHOLD.
                                 // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_H_ANISO;  // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See the Used Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double SPACE_ERR_STOP = 5.0;       // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit. This is mainly to prevent h-adaptivity to go on forever.

// Problem parameters
const double RE = 200.0;             // Reynolds number.
const double VEL_INLET = 1.0;        // Inlet velocity (reached after STARTUP_TIME).
const double STARTUP_TIME = 1.0;     // During this time, inlet velocity increases gradually
                                     // from 0 to VEL_INLET, then it stays constant.
const double TAU = 0.01;             // Time step.
const double T_FINAL = 30000.0;      // Time interval length.

// Newton's method
const double NEWTON_TOL_COARSE = 0.01;     // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.05;       // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 20;            // Maximum allowed number of Newton iterations.

// Geometry
const double H = 5;                  // Domain height (necessary to define the parabolic
                                     // velocity profile at inlet)

// Boundary markers.
int bdy_bottom = 1;
int bdy_right  = 2;
int bdy_top = 3;
int bdy_left = 4;
int bdy_obstacle = 5;

// Current time (defined as global since needed in weak forms)
double TIME = 0;

// Boundary condition types for x-velocity
BCType xvel_bc_type(int marker) {
  if (marker == bdy_right) return BC_NONE;
  else return BC_ESSENTIAL;
}

// Boundary condition values for x-velocity
scalar essential_bc_values_xvel(int ess_bdy_marker, double x, double y) {
  if (ess_bdy_marker == bdy_left) {
    // time-dependent inlet velocity (parabolic profile)
    double val_y = VEL_INLET * y*(H-y) / (H/2.)/(H/2.); //parabolic profile with peak VEL_INLET at y = H/2
    if (TIME <= STARTUP_TIME) return val_y * TIME/STARTUP_TIME;
    else return val_y;
  }
  else return 0;
}

// Essential (Dirichlet) boundary condition values for y-velocity.
scalar essential_bc_values_yvel(int ess_bdy_marker, double x, double y) 
{
  return 0;
}

// Boundary condition types for y-velocity
BCType yvel_bc_type(int marker) {
  if (marker == bdy_right) return BC_NONE;
  else return BC_ESSENTIAL;
}

BCType p_bc_type(int marker)
  { return BC_NONE; }

// Weak forms
#include "forms.cpp"

void mag(int n, scalar* a, scalar* dadx, scalar* dady,
                scalar* b, scalar* dbdx, scalar* dbdy,
                scalar* out, scalar* outdx, scalar* outdy)
{
  for (int i = 0; i < n; i++)
  {
    out[i] = sqrt(sqr(a[i]) + sqr(b[i]));
    outdx[i] = (0.5 / out[i]) * (2.0 * a[i] * dadx[i] + 2.0 * b[i] * dbdx[i]);
    outdy[i] = (0.5 / out[i]) * (2.0 * a[i] * dady[i] + 2.0 * b[i] * dbdy[i]);
  }
}

int main(int argc, char* argv[])
{
  // Load the mesh file.
  Mesh basemesh, mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &basemesh);  // Master mesh.

  // Perform initial mesh refinements.
  for (int i=0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  basemesh.refine_towards_boundary(bdy_obstacle, INIT_REF_NUM_BDY, false); // 'true' stands for anisotropic refinements,
  basemesh.refine_towards_boundary(bdy_top, INIT_REF_NUM_BDY, true);       // 'false' for isotropic.
  basemesh.refine_towards_boundary(bdy_bottom, INIT_REF_NUM_BDY, true);
  mesh.copy(&basemesh);

  // Create spaces with default shapesets. 
  H1Space xvel_space(&mesh, xvel_bc_type, essential_bc_values_xvel, P_INIT_VEL);
  H1Space yvel_space(&mesh, yvel_bc_type, essential_bc_values_yvel, P_INIT_VEL);
#ifdef PRESSURE_IN_L2
  L2Space p_space(&mesh, P_INIT_PRESSURE);
#else
  H1Space p_space(&mesh, p_bc_type, NULL, P_INIT_PRESSURE);
#endif

  // Define projection norms.
  int vel_proj_norm = 1;
#ifdef PRESSURE_IN_L2
  int p_proj_norm = 0;
#else
  int p_proj_norm = 1;
#endif

  // Solutions for the Newton's iteration and time stepping.
  Solution xvel_fine, yvel_fine, p_fine;
  Solution xvel_prev_time, yvel_prev_time, p_prev_time;
  Solution xvel_prev_newton, yvel_prev_newton, p_prev_newton;

  // Define initial conditions on the coarse mesh.
  xvel_prev_time.set_zero(&mesh);
  yvel_prev_time.set_zero(&mesh);
  p_prev_time.set_zero(&mesh);

  // Initialize the weak formulation.
  WeakForm wf(3);
  wf.add_matrix_form(0, 0, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
  wf.add_matrix_form(0, 0, callback(newton_bilinear_form_unsym_0_0), H2D_UNSYM, H2D_ANY, 
                Tuple<MeshFunction*>(&xvel_prev_newton, &yvel_prev_newton));
  wf.add_matrix_form(0, 1, callback(newton_bilinear_form_unsym_0_1), H2D_UNSYM, H2D_ANY, 
                &xvel_prev_newton);
  wf.add_matrix_form(0, 2, callback(bilinear_form_unsym_0_2), H2D_ANTISYM);
  wf.add_matrix_form(1, 0, callback(newton_bilinear_form_unsym_1_0), H2D_UNSYM, H2D_ANY, 
                &yvel_prev_newton);
  wf.add_matrix_form(1, 1, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
  wf.add_matrix_form(1, 1, callback(newton_bilinear_form_unsym_1_1), H2D_UNSYM, H2D_ANY, 
                Tuple<MeshFunction*>(&xvel_prev_newton, &yvel_prev_newton));
  wf.add_matrix_form(1, 2, callback(bilinear_form_unsym_1_2), H2D_ANTISYM);
  wf.add_vector_form(0, callback(newton_F_0), H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time, 
							 &xvel_prev_newton, &yvel_prev_newton, &p_prev_newton));
  wf.add_vector_form(1, callback(newton_F_1), H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time, 
							 &xvel_prev_newton, &yvel_prev_newton, &p_prev_newton));
  wf.add_vector_form(2, callback(newton_F_2), H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_newton, &yvel_prev_newton));

  // Initialize nonlinear system.
  NonlinSystem nls(&wf, Tuple<Space*>(&xvel_space, &yvel_space, &p_space));

  // Update time-dependent essential BC.
  info("Updating time-dependent essential BC.");
  TIME += TAU;
  nls.update_essential_bc_values();

  // Project initial conditions on FE spaces to obtain initial coefficient 
  // vector for the Newton's method.
  info("Projecting initial conditions to obtain initial vector for the Newton'w method.");
  nls.project_global(Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time, &p_prev_time),
                     Tuple<Solution*>(&xvel_prev_newton, &yvel_prev_newton, &p_prev_newton),
                     Tuple<int>(vel_proj_norm, vel_proj_norm, p_proj_norm));  

  // View the initial condition.
  VectorView vview("Velocity initial condition", 0, 0, 750, 240);
  ScalarView pview("Pressure initial condition", 0, 290, 750, 240);
  vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  //pview.set_min_max_range(-0.9, 1.0);
  pview.fix_scale_width(80);
  pview.show_mesh(true);
  vview.show(&xvel_prev_time, &yvel_prev_time, H2D_EPS_LOW);
  pview.show(&p_prev_time);

  // Solve on the coarse meshes.
  // Newton's loop on the coarse mesh.
  bool verbose = true; // Default is false.
  info("Solving on coarse meshes.");
  if (!nls.solve_newton(Tuple<Solution*>(&xvel_prev_newton, &yvel_prev_newton, &p_prev_newton), 
                        NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
    error("Newton's method did not converge.");

  // Store the results on coarse meshes.
  Solution xvel_coarse, yvel_coarse, p_coarse;
  xvel_coarse.copy(&xvel_prev_newton);
  yvel_coarse.copy(&yvel_prev_newton);
  p_coarse.copy(&p_prev_newton);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Time-stepping loop:
  int num_time_steps = (int)(T_FINAL / TAU + 0.5);
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      xvel_space.set_uniform_order(P_INIT_VEL);
      yvel_space.set_uniform_order(P_INIT_VEL);
      p_space.set_uniform_order(P_INIT_PRESSURE);

      // Project fine mesh solution on the globally derefined meshes.
      info("---- Time step %d:", ts);
      if (SOLVE_ON_COARSE_MESH) 
        info("Projecting fine mesh solutions to obtain initial vector on globally derefined meshes.");
      else 
        info("Projecting fine mesh solutions on globally derefined meshes for error calculation.");
      nls.project_global(Tuple<MeshFunction*>(&xvel_prev_newton, &yvel_prev_newton, &p_prev_newton),
                         Tuple<Solution*>(&xvel_prev_newton, &yvel_prev_newton, &p_prev_newton),
                         Tuple<int>(vel_proj_norm, vel_proj_norm, p_proj_norm));

      if (SOLVE_ON_COARSE_MESH) {
        // Newton's loop on the globally derefined meshes.
        info("Solving on globally derefined meshes.", ts);
        if (!nls.solve_newton(Tuple<Solution*>(&xvel_prev_newton, &yvel_prev_newton, &p_prev_newton), 
                              NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
          error("Newton's method did not converge.");
      }

      // Store the result on coarse meshes..
      xvel_coarse.copy(&xvel_prev_newton);
      yvel_coarse.copy(&yvel_prev_newton);
      p_coarse.copy(&p_prev_newton);
    }

    // Adaptivity loop (in space):
    bool done = false;
    double space_err_est;
    int as = 1;

    do
    {
      info("---- Time step %d, adaptivity step %d:", ts, as);

      // Initialize reference nonlinear system.
      int p_increment = 0;
      RefSystem rnls(&nls, p_increment);

      // Set initial condition for the Newton's method on the fine meshes.
      if (as == 1) {
        info("Projecting coarse mesh solutions to obtain initial vector on new fine meshes.");
        rnls.project_global(Tuple<MeshFunction*>(&xvel_coarse, &yvel_coarse, &p_coarse),
                            Tuple<Solution*>(&xvel_prev_newton, &yvel_prev_newton, &p_prev_newton),
                            Tuple<int>(vel_proj_norm, vel_proj_norm, p_proj_norm));
      }
      else {
        info("Projecting fine mesh solution to obtain initial vector on new fine mesh.");
        rnls.project_global(Tuple<MeshFunction*>(&xvel_fine, &yvel_fine, &p_fine),
                            Tuple<Solution*>(&xvel_prev_newton, &yvel_prev_newton, &p_prev_newton),
                            Tuple<int>(vel_proj_norm, vel_proj_norm, p_proj_norm));
      }

      // Newton's method on fine meshes
      info("Solving on fine meshes.");
      if (!rnls.solve_newton(Tuple<Solution*>(&xvel_prev_newton, &yvel_prev_newton, &p_prev_newton),
                             NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose))
        error("Newton's method did not converge.");

      // Store the result on fine meshes.
      xvel_fine.copy(&xvel_prev_newton);
      yvel_fine.copy(&yvel_prev_newton);
      p_fine.copy(&p_prev_newton);

      /*
      // Measure error in velocity magnitude.
      DXDYFilter crs_mag(mag, &xvel_crs, &yvel_crs);
      DXDYFilter fine_mag(mag, &xvel_prev_newton, &yvel_prev_newton);
      double space_err = 100 * h1_error(&crs_mag, &fine_mag);
      info("Velocity rel error est %g%%", space_err);
      */

      H1Adapt hp(&nls);
      hp.set_solutions(Tuple<Solution*>(&xvel_coarse, &yvel_coarse, &p_coarse),
                       Tuple<Solution*>(&xvel_fine, &yvel_fine, &p_fine));
      // Error calculated using velocity components only.
      hp.set_biform(0, 0, callback(h1_form));
      hp.set_biform(1, 1, callback(h1_form));
      double space_err_est = hp.calc_error() * 100;
      info("ndof_coarse: %d, ndof_fine: %d, space_err_est: %g%%", 
        nls.get_num_dofs(), rnls.get_num_dofs(), space_err_est);
 
      // If space_err_est too large, adapt the mesh.
      if (space_err_est < SPACE_ERR_STOP) done = true;
      else {
        info("Adapting coarse meshes.");
        done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
	printf("aaa\n");
        if (nls.get_num_dofs() >= NDOF_STOP) {
          done = true;
          break;
        }
	printf("bbb\n");

        // Project the fine mesh solution on the new coarse mesh.
        if (SOLVE_ON_COARSE_MESH) 
          info("Projecting fine mesh solutions to obtain initial vector on new coarse meshes.");
        else 
          info("Projecting fine mesh solutions on coarse meshes for error calculation.");
        nls.project_global(Tuple<MeshFunction*>(&xvel_fine, &yvel_fine, &p_fine),
                           Tuple<Solution*>(&xvel_prev_newton, &yvel_prev_newton, &p_prev_newton),
                           Tuple<int>(vel_proj_norm, vel_proj_norm, p_proj_norm));

        if (SOLVE_ON_COARSE_MESH) {
          // Newton's loop on the coarse meshes.
          info("---- Time step %d, adaptivity step %d, solving on new coarse mesh.", ts, as);
          if (!nls.solve_newton(Tuple<Solution*>(&xvel_prev_newton, &yvel_prev_newton, &p_prev_newton), 
                                NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
            error("Newton's method did not converge.");
        }

        // Store the result in sln_coarse.
        xvel_coarse.copy(&xvel_prev_newton);
        yvel_coarse.copy(&yvel_prev_newton);
        p_coarse.copy(&p_prev_newton);
      
        as++;
      }
    }
    while (!done);

    // Save solutions for the next time step.
    xvel_prev_time.copy(&xvel_fine);
    yvel_prev_time.copy(&yvel_fine);
    p_prev_time.copy(&p_fine);

    // Visualize the solution and mesh.
    char title[100];
    sprintf(title, "Velocity, time %g s", TIME);
    vview.set_title(title);
    vview.show(&xvel_coarse, &yvel_coarse, H2D_EPS_LOW);
    sprintf(title, "Pressure, time %g s", TIME);
    pview.set_title(title);
    pview.show(&p_coarse);

    // Update global time.
    TIME += TAU;

    // Update time dependent essential BC.
    if (TIME <= STARTUP_TIME) {
      info("Updating time-dependent essential BC.");
      nls.update_essential_bc_values();
    }
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
