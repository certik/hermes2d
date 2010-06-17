#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

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

//#define STOKES                     // If this is defined, Stokes problem is solved, otherwise N-S.
#define PRESSURE_IN_L2               // If this is defined, the pressure is approximated using
                                     // discontinuous L2 elements (making the velocity discreetely
                                     // divergence-free, more accurate than using a continuous
                                     // pressure approximation). Otherwise the standard continuous
                                     // elements are used. The results are striking - check the
                                     // tutorial for comparisons.
const int INIT_REF_NUM = 1;          // Number of initial uniform mesh refinements. 
const int INIT_BDY_REF_NUM = 1;      // Number of initial mesh refinements towards boundary. 
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
const double TAU = 0.1;             // Time step.
const double T_FINAL = 30000.0;      // Time interval length.
const double NEWTON_TOL = 1e-5;      // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 10;      // Maximum allowed number of Newton iterations.

// Boundary markers.
int bdy_inner = 1;
int bdy_outer  = 2;

// Current time (used in weak forms).
double TIME = 0;

// Boundary condition types for velocity.
BCType vel_bc_type(int marker) {
  return BC_ESSENTIAL;
}

// Boundary condition types for pressure.
BCType p_bc_type(int marker)
  { return BC_NONE; }

// Essential (Dirichlet) boundary condition values for x-velocity.
scalar essential_bc_value_xvel(int ess_bdy_marker, double x, double y) {
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
scalar essential_bc_value_yvel(int ess_bdy_marker, double x, double y) {
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

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Initial mesh refinements.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(bdy_inner, INIT_BDY_REF_NUM, false); // INIT_BDY_REF_NUM is the number of levels,
  //mesh.refine_towards_boundary(bdy_outer, INIT_BDY_REF_NUM, true); // 'true' stands for anisotropic refinements.

  // Spaces for velocity components and pressure.
  H1Space xvel_space(&mesh, vel_bc_type, essential_bc_value_xvel, P_INIT_VEL);
  H1Space yvel_space(&mesh, vel_bc_type, essential_bc_value_yvel, P_INIT_VEL);
#ifdef PRESSURE_IN_L2
  L2Space p_space(&mesh, P_INIT_VEL);
#else
  H1Space p_space(&mesh, p_bc_type, NULL, P_INIT_PRESSURE);
#endif

  // Solutions for the Newton's iteration and time stepping.
  Solution xvel_prev_time, yvel_prev_time, xvel_prev_newton, yvel_prev_newton, p_prev;
  xvel_prev_time.set_zero(&mesh);
  yvel_prev_time.set_zero(&mesh);
  xvel_prev_newton.set_zero(&mesh);
  yvel_prev_newton.set_zero(&mesh);
  p_prev.set_zero(&mesh);

  // Initialize weak formulation.
  WeakForm wf(3);
  if (NEWTON) {
    wf.add_biform(0, 0, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_biform(0, 0, callback(newton_bilinear_form_unsym_0_0), 
                  H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_newton, &yvel_prev_newton));
    wf.add_biform(0, 1, callback(newton_bilinear_form_unsym_0_1), 
                  H2D_UNSYM, H2D_ANY, &xvel_prev_newton);
    wf.add_biform(0, 2, callback(bilinear_form_unsym_0_2), H2D_ANTISYM);
    wf.add_biform(1, 0, callback(newton_bilinear_form_unsym_1_0), 
                  H2D_UNSYM, H2D_ANY, &yvel_prev_newton);
    wf.add_biform(1, 1, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_biform(1, 1, callback(newton_bilinear_form_unsym_1_1), 
                  H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_newton, &yvel_prev_newton));
    wf.add_biform(1, 2, callback(bilinear_form_unsym_1_2), H2D_ANTISYM);
    wf.add_liform(0, callback(newton_F_0), H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_time, 
									 &yvel_prev_time, &xvel_prev_newton, &yvel_prev_newton, &p_prev));
    wf.add_liform(1, callback(newton_F_1), H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_time, 
									 &yvel_prev_time, &xvel_prev_newton, &yvel_prev_newton, &p_prev));
    wf.add_liform(2, callback(newton_F_2), H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_newton, &yvel_prev_newton));
  }
  else {
    wf.add_biform(0, 0, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_biform(0, 0, callback(simple_bilinear_form_unsym_0_0_1_1), 
                  H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time));
    wf.add_biform(1, 1, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_biform(1, 1, callback(simple_bilinear_form_unsym_0_0_1_1), 
                  H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time));
    wf.add_biform(0, 2, callback(bilinear_form_unsym_0_2), H2D_ANTISYM);
    wf.add_biform(1, 2, callback(bilinear_form_unsym_1_2), H2D_ANTISYM);
    wf.add_liform(0, callback(simple_linear_form), H2D_ANY, &xvel_prev_time);
    wf.add_liform(1, callback(simple_linear_form), H2D_ANY, &yvel_prev_time);
  }

  // Initialize views.
  VectorView vview("velocity [m/s]", 0, 0, 500, 400);
  ScalarView pview("pressure [Pa]", 510, 0, 500, 400);
  //vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  //pview.set_min_max_range(-0.9, 1.0);
  pview.fix_scale_width(80);
  pview.show_mesh(true);

  // Matrix solver.
  UmfpackSolver umfpack;

  // Initialize linear system.
  /*
  LinSystem ls(&wf, &umfpack, 3, &xvel_space, &yvel_space, &p_space);
  if (!NEWTON) printf("LinSystem ndof = %d\n", ls.get_num_dofs());
  if (!NEWTON) ls.project_global(&xvel_prev_newton, &yvel_prev_newton, &p_prev, 
                                 &xvel_prev_newton, &yvel_prev_newton, &p_prev);
  */

  // Initialize nonlinear system.
  NonlinSystem nls(&wf, &umfpack, Tuple<Space*>(&xvel_space, &yvel_space, &p_space));

  // View FE basis functions.
  //BaseView view1;
  //view1.show(&p_space);
  //View::wait();

  if (NEWTON) printf("NonlinSystem ndof = %d\n", nls.get_num_dofs());
  if(nls.get_solution_vector() == NULL) error("Vec is NULL.");
  if (NEWTON) nls.project_global(Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time, &p_prev), 
                                 Tuple<Solution*>(&xvel_prev_time, &yvel_prev_time, &p_prev));

  // Setting initial conditions for Newton's iteration.
  xvel_prev_newton.copy(&xvel_prev_time);
  yvel_prev_newton.copy(&yvel_prev_time);

  // Show initial condition for velocity. 
  vview.show(&xvel_prev_time, &yvel_prev_time, H2D_EPS_LOW);
  vview.wait(H2DV_WAIT_KEYPRESS);

  // Time-stepping loop:
  char title[100];
  int num_time_steps = T_FINAL / TAU;
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    TIME += TAU;

    info("---- Time step %d, time = %g:", ts, TIME);

    // This is needed to update time-dependent boundary conditions.
    update_essential_bc_values(Tuple<Space*>(&xvel_space, &yvel_space, &p_space));

    if (NEWTON) {
      // Newton's method.
      info("Performing Newton's method.");
      if (!nls.solve_newton(&xvel_prev_newton, &yvel_prev_newton, &p_prev, NEWTON_TOL, NEWTON_MAX_ITER)) 
        error("Newton's method did not converge.");

      // Show the solution at the end of time step.
      sprintf(title, "Velocity, time %g", TIME);
      vview.set_title(title);
      vview.show(&xvel_prev_newton, &yvel_prev_newton, H2D_EPS_LOW);
      sprintf(title, "Pressure, time %g", TIME);
      pview.set_title(title);
      pview.show(&p_prev);

      // Calculate an estimate of the temporal change of the x-velocity.
      H1Adapt hp(&xvel_space);
      hp.set_solutions(&xvel_prev_time, &xvel_prev_newton);
      double err_est = hp.calc_error(H2D_TOTAL_ERROR_ABS | H2D_ELEMENT_ERROR_ABS) / TAU;
      info("x_vel temporal change: %g", err_est);

      // Copy the result of the Newton's iteration into the
      // previous time level solutions.
      xvel_prev_time.copy(&xvel_prev_newton);
      yvel_prev_time.copy(&yvel_prev_newton);
    }
    else {
      // Assemble and solve.
      Solution xvel_sln, yvel_sln, p_sln;
      info("Assembling and solving linear problem.");
      /*
      ls.assemble();
      ls.solve(3, &xvel_sln, &yvel_sln, &p_sln);
      */

      // Show the solution at the end of time step.
      sprintf(title, "Velocity, time %g", TIME);
      vview.set_title(title);
      vview.show(&xvel_sln, &yvel_sln, H2D_EPS_LOW);
      sprintf(title, "Pressure, time %g", TIME);
      pview.set_title(title);
      pview.show(&p_sln);

      // Calculate an estimate of the temporal change of the x-velocity.
      H1Adapt hp(&xvel_space);
      hp.set_solutions(&xvel_prev_time, &xvel_sln);
      double err_est = hp.calc_error(H2D_TOTAL_ERROR_ABS | H2D_ELEMENT_ERROR_ABS) / TAU;
      info("x_vel temporal change: %g", err_est);
  
      // This copy destroys xvel_sln and yvel_sln
      // which are no longer needed.
      xvel_prev_time = xvel_sln;
      yvel_prev_time = yvel_sln;
    }
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
