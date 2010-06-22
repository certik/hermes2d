#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "forms.h"

// Flow in between two circles, inner circle is rotating with surface 
// velocity VEL. The time-dependent laminar incompressible Navier-Stokes equations
// are discretized in time via the implicit Euler method. If NEWTON == true,
// the Newton's method is used to solve the nonlinear problem at each time
// step. If NEWTON == false, the convective term is only linearized using the
// velocities from the previous time step. Obviously the latter approach is wrong,
// but people do this frequently because it is faster and simpler to implement.
// Therefore we include this case for comparison purposes. We also show how
// to use discontinuous ($L^2$) elements for pressure and thus make the
// velocity discretely divergence-free. Comparison to approximating the
// pressure with the standard (continuous) Taylor-Hood elements is enabled.
// The Reynolds number Re = 200 which is very low. You can increase it but 
// then you will need to make the mesh finer, and the computation will take 
// more time.
//
// PDE: incompressible Navier-Stokes equations in the form
//     \partial v / \partial t - \Delta v / Re + (v \cdot \nabla) v + \nabla p = 0,
//     div v = 0.
//
// BC: tangential velocity V on Gamma_1 (inner circle),
//     zero velocity on Gamma_2 (outer circle).
//
// Geometry: Area in between two concentric circles with radiuses r1 and r2,
//           r1 < r2. These radiuses can be changed in the file "domain.mesh".
//
// The following parameters can be changed:

const int INIT_REF_NUM = 2;                // Number of initial uniform mesh refinements. 
const int INIT_BDY_REF_NUM_INNER = 2;      // Number of initial mesh refinements towards boundary. 
const int INIT_BDY_REF_NUM_OUTER = 2;      // Number of initial mesh refinements towards boundary. 

//#define STOKES                     // If this is defined, Stokes problem is solved, otherwise N-S.
#define PRESSURE_IN_L2               // If this is defined, the pressure is approximated using
                                     // discontinuous L2 elements (making the velocity discreetely
                                     // divergence-free, more accurate than using a continuous
                                     // pressure approximation). Otherwise the standard continuous
                                     // elements are used. The results are striking - check the
                                     // tutorial for comparisons.
const bool NEWTON = true;            // If NEWTON == true then the Newton's iteration is performed.
                                     // in every time step. Otherwise the convective term is linearized
                                     // using the velocities from the previous time step.
const int P_INIT_VEL = 2;            // Initial polynomial degree for velocity components.
const int P_INIT_PRESSURE = 1;       // Initial polynomial degree for pressure.
                                     // Note: P_INIT_VEL should always be greater than
                                     // P_INIT_PRESSURE because of the inf-sup condition.
const double RE = 200.0;             // Reynolds number.
const double VEL = 0.1;              // Surface velocity of inner circle.
const double STARTUP_TIME = 1.0;     // During this time, surface velocity of the inner circle increases 
                                     // gradually from 0 to VEL, then it stays constant.
const double TAU = 0.1;              // Time step.
const double T_FINAL = 30000.0;      // Time interval length.
const double NEWTON_TOL = 1e-5;      // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 10;      // Maximum allowed number of Newton iterations.

// Boundary markers.
int bdy_inner = 1;
int bdy_outer = 2;

// Current time (used in weak forms).
double TIME = 0;

// Boundary condition types for x-velocity.
BCType xvel_bc_type(int marker) {
  return BC_ESSENTIAL;
}

// Boundary condition types for y-velocity.
BCType yvel_bc_type(int marker) {
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values for x-velocity.
scalar essential_bc_values_xvel(int ess_bdy_marker, double x, double y) {
  if (ess_bdy_marker == bdy_inner) {
    // Time-dependent surface velocity of inner circle.
    double velocity;
    if (TIME <= STARTUP_TIME) velocity = VEL * TIME/STARTUP_TIME;
    else velocity = VEL;
    double alpha = atan2(x, y);
    double xvel = velocity*cos(alpha);
    //printf("%g %g xvel = %g\n", x, y, xvel);
    return xvel; 
  }
  else return 0;
}

// Essential (Dirichlet) boundary condition values for y-velocity.
scalar essential_bc_values_yvel(int ess_bdy_marker, double x, double y) {
  if (ess_bdy_marker == bdy_inner) {
    // Time-dependent surface velocity of inner circle.
    double velocity;
    if (TIME <= STARTUP_TIME) velocity = VEL * TIME/STARTUP_TIME;
    else velocity = VEL;
    double alpha = atan2(x, y);
    double yvel = -velocity*sin(alpha);
    //printf("%g %g yvel = %g\n", x, y, yvel);
    return yvel; 
  }
  else return 0;
}

// Weak forms.
#include "forms.cpp"

// Custom function to calculate drag coefficient.
double integrate_over_wall(MeshFunction* meshfn, int marker)
{
  Quad2D* quad = &g_quad_2d_std;
  meshfn->set_quad_2d(quad);

  double integral = 0.0;
  Element* e;
  Mesh* mesh = meshfn->get_mesh();

  for_all_active_elements(e, mesh)
  {
    for(int edge = 0; edge < e->nvert; edge++)
    {
      if ((e->en[edge]->bnd) && (e->en[edge]->marker == marker))
      {
        update_limit_table(e->get_mode());
        RefMap* ru = meshfn->get_refmap();

        meshfn->set_active_element(e);
        int eo = quad->get_edge_points(edge);
        meshfn->set_quad_order(eo, H2D_FN_VAL);
        scalar *uval = meshfn->get_fn_values();
        double3* pt = quad->get_points(eo);
        double3* tan = ru->get_tangent(edge);
        for (int i = 0; i < quad->get_num_points(eo); i++)
          integral += pt[i][2] * uval[i] * tan[i][2];
      }
    }
  }
  return integral * 0.5;
}

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain-excentric.mesh", &mesh);
  //mloader.load("domain-concentric.mesh", &mesh);

  // Initial mesh refinements.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(bdy_inner, INIT_BDY_REF_NUM_INNER, false);  // true for anisotropic refinements
  mesh.refine_towards_boundary(bdy_outer, INIT_BDY_REF_NUM_OUTER, false);  // false for isotropic refinements

  // View the mesh.
  //MeshView mv("Mesh", 0, 0, 500, 500);
  //mv.show(&mesh);
  //View::wait();

  // Create spaces with default shapesets. 
  H1Space xvel_space(&mesh, xvel_bc_type, essential_bc_values_xvel, P_INIT_VEL);
  H1Space yvel_space(&mesh, yvel_bc_type, essential_bc_values_yvel, P_INIT_VEL);
#ifdef PRESSURE_IN_L2
  L2Space p_space(&mesh, P_INIT_PRESSURE);
#else
  H1Space p_space(&mesh, NULL, NULL, P_INIT_PRESSURE);
#endif

  // Define projection norms.
  int vel_proj_norm = 1;
#ifdef PRESSURE_IN_L2
  int p_proj_norm = 0;
#else
  int p_proj_norm = 1;
#endif

  // Solutions for the Newton's iteration and time stepping.
  Solution xvel_prev_time, yvel_prev_time, p_prev_time;
  Solution xvel_prev_newton, yvel_prev_newton, p_prev_newton;

  // Initialize weak formulation.
  WeakForm wf(3);
  if (NEWTON) {
    wf.add_matrix_form(0, 0, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_matrix_form(0, 0, callback(newton_bilinear_form_unsym_0_0), 
                  H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_newton, &yvel_prev_newton));
    wf.add_matrix_form(0, 1, callback(newton_bilinear_form_unsym_0_1), 
                  H2D_UNSYM, H2D_ANY, &xvel_prev_newton);
    wf.add_matrix_form(0, 2, callback(bilinear_form_unsym_0_2), H2D_ANTISYM);
    wf.add_matrix_form(1, 0, callback(newton_bilinear_form_unsym_1_0), 
                  H2D_UNSYM, H2D_ANY, &yvel_prev_newton);
    wf.add_matrix_form(1, 1, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_matrix_form(1, 1, callback(newton_bilinear_form_unsym_1_1), 
                  H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_newton, &yvel_prev_newton));
    wf.add_matrix_form(1, 2, callback(bilinear_form_unsym_1_2), H2D_ANTISYM);
    wf.add_vector_form(0, callback(newton_F_0), H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_time, 
		  &yvel_prev_time, &xvel_prev_newton, &yvel_prev_newton, &p_prev_newton));
    wf.add_vector_form(1, callback(newton_F_1), H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_time, 
		  &yvel_prev_time, &xvel_prev_newton, &yvel_prev_newton, &p_prev_newton));
    wf.add_vector_form(2, callback(newton_F_2), H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_newton, 
                  &yvel_prev_newton));
  }
  else {
    wf.add_matrix_form(0, 0, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_matrix_form(0, 0, callback(simple_bilinear_form_unsym_0_0_1_1), 
                  H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time));
    wf.add_matrix_form(1, 1, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_matrix_form(1, 1, callback(simple_bilinear_form_unsym_0_0_1_1), 
                  H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time));
    wf.add_matrix_form(0, 2, callback(bilinear_form_unsym_0_2), H2D_ANTISYM);
    wf.add_matrix_form(1, 2, callback(bilinear_form_unsym_1_2), H2D_ANTISYM);
    wf.add_vector_form(0, callback(simple_linear_form), H2D_ANY, &xvel_prev_time);
    wf.add_vector_form(1, callback(simple_linear_form), H2D_ANY, &yvel_prev_time);
  }

  // Initialize views.
  VectorView vview("velocity [m/s]", 0, 0, 600, 500);
  ScalarView pview("pressure [Pa]", 610, 0, 600, 500);
  //vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  //pview.set_min_max_range(-0.9, 1.0);
  pview.fix_scale_width(80);
  pview.show_mesh(true);

  // Initialize linear system.
  LinSystem ls(&wf, Tuple<Space*>(&xvel_space, &yvel_space, &p_space));

  // Initialize nonlinear system.
  NonlinSystem nls(&wf, Tuple<Space*>(&xvel_space, &yvel_space, &p_space));

  // Set initial conditions.
  info("Setting initial conditions.");
  xvel_prev_time.set_zero(&mesh);
  yvel_prev_time.set_zero(&mesh);
  p_prev_time.set_zero(&mesh);

  // Update time-dependent essential BC.
  info("Updating time-dependent essential BC.");
  TIME += TAU;
  if (NEWTON) nls.update_essential_bc_values();
  else ls.update_essential_bc_values();

  // Project initial conditions on the FE spaces
  // to obtain initial coefficient vector for the Newton's method.
  info("Projecting initial conditions to obtain initial vector for the Newton'w method.");
  nls.project_global(Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time, &p_prev_time),
                     Tuple<Solution*>(&xvel_prev_newton, &yvel_prev_newton, &p_prev_newton),
                     Tuple<int>(vel_proj_norm, vel_proj_norm, p_proj_norm));  

  // Time-stepping loop:
  char title[100];
  int num_time_steps = T_FINAL / TAU;
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    TIME += TAU;
    info("---- Time step %d, time = %g:", ts, TIME);

    if (NEWTON) {
      // Newton's method.
      info("Performing Newton's method.");
      bool verbose = true; // Default is false.
      if (!nls.solve_newton(Tuple<Solution*>(&xvel_prev_newton, &yvel_prev_newton, &p_prev_newton), 
                            NEWTON_TOL, NEWTON_MAX_ITER, verbose)) {
        error("Newton's method did not converge.");
      }
    }
    else {
      // Assemble and solve.
      info("Assembling and solving linear problem.");
      ls.assemble();
      ls.solve(Tuple<Solution*>(&xvel_prev_newton, &yvel_prev_newton, &p_prev_newton));
    }

    // Calculate an estimate of the temporal change of the velocity.
    H1Adapt hp(&nls);
    hp.set_solutions(Tuple<Solution*>(&xvel_prev_time, &yvel_prev_time, &p_prev_time), 
                     Tuple<Solution*>(&xvel_prev_newton, &yvel_prev_newton, &p_prev_newton));
    hp.set_biform(0, 0, callback(h1_form));
    hp.set_biform(1, 1, callback(h1_form));
    double err_est = hp.calc_error(H2D_TOTAL_ERROR_ABS | H2D_ELEMENT_ERROR_ABS) / TAU;
    info("x_vel temporal change: %g", err_est);
    // Show the solution at the end of time step.
    sprintf(title, "Velocity, time %g", TIME);
    vview.set_title(title);
    vview.show(&xvel_prev_newton, &yvel_prev_newton, H2D_EPS_LOW);
    sprintf(title, "Pressure, time %g", TIME);
    pview.set_title(title);
    pview.show(&p_prev_newton);

    // Calculate pressure integral over inner circle.
    // (To be replaced with drag coefficient calculation later.)
    double val = integrate_over_wall(&p_prev_newton, bdy_inner);    
    printf("Pressure integral: %g\n", val);

    // Copy the result of the Newton's iteration into the
    // previous time level solutions.
    xvel_prev_time.copy(&xvel_prev_newton);
    yvel_prev_time.copy(&yvel_prev_newton);
    p_prev_time.copy(&p_prev_newton);

    // Update global time.
    TIME += TAU;

    // Update time dependent essential BC.
    if (TIME <= STARTUP_TIME) {
      info("Updating time-dependent essential BC.");
      if (NEWTON) nls.update_essential_bc_values();
      else ls.update_essential_bc_values();
    }
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
