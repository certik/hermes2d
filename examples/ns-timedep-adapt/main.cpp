#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

// The time-dependent laminar incompressible Navier-Stokes equations are
// discretized in time via the implicit Euler method. If NEWTON == true,
// the Newton's method is used to solve the nonlinear problem at each time
// step. If NEWTON == false, the convective term is only linearized using the
// velocities from the previous time step. Obviously the latter approach is wrong,
// but people do this frequently because it is faster and simpler to implement.
// Therefore we include this case for comparison purposes. We also show how
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
//
// The following parameters can be changed:
//

#define PRESSURE_IN_L2           // If this is defined, the pressure is approximated using
                                 // discontinuous L2 elements (making the velocity discreetely
                                 // divergence-free, more accurate than using a continuous
                                 // pressure approximation). Otherwise the standard continuous
                                 // elements are used. The results are striking - check the
                                 // tutorial for comparisons.
const bool NEWTON = true;        // If NEWTON == true then the Newton's iteration is performed
                                 // in every time step. Otherwise the convective term is linearized
                                 // using the velocities from the previous time step
const int P_INIT_VEL = 2;        // Initial polynomial degree for velocity components
const int P_INIT_PRESSURE = 1;   // Initial polynomial degree for pressure
                                 // Note: P_INIT_VEL should always be greater than
                                 // P_INIT_PRESSURE because of the inf-sup condition
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
const int ADAPT_TYPE = 1;        // Type of automatic adaptivity:
                                 // ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                 // ADAPT_TYPE = 1 ... adaptive h-FEM,
                                 // ADAPT_TYPE = 2 ... adaptive p-FEM.
const bool ISO_ONLY = false;     // Isotropic refinement flag (concerns quadrilateral elements only).
                                 // ISO_ONLY = false ... anisotropic refinement of quad elements
                                 // is allowed (default),
                                 // ISO_ONLY = true ... only isotropic refinements of quad elements
                                 // are allowed.
const int MESH_REGULARITY = -1;  // Maximum allowed level of hanging nodes:
                                 // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                 // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                 // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                 // Note that regular meshes are not supported, this is due to
                                 // their notoriously bad performance.
const double CONV_EXP = 1.0;     // Default value is 1.0. This parameter influences the selection of
                                 // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double SPACE_TOL = 1.0;    // Stopping criterion for adaptivity (rel. error tolerance between the
                                 // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;     // Adaptivity process stops when the number of degrees of freedom grows over
                                 // this limit. This is mainly to prevent h-adaptivity to go on forever.

// Problem parameters
const double RE = 200.0;             // Reynolds number
const double VEL_INLET = 1.0;        // Inlet velocity (reached after STARTUP_TIME)
const double STARTUP_TIME = 1.0;     // During this time, inlet velocity increases gradually
                                     // from 0 to VEL_INLET, then it stays constant
const double TAU = 0.1;              // time step
const double T_FINAL = 30000.0;      // Time interval length
const double NEWTON_TOL = 1e-3;      // Stopping criterion for the Newton's method
const int NEWTON_MAX_ITER = 10;      // Maximum allowed number of Newton iterations
const double H = 5;                  // Domain height (necessary to define the parabolic
                                     // velocity profile at inlet)

// Boundary markers in the mesh file
int marker_bottom = 1;
int marker_right  = 2;
int marker_top = 3;
int marker_left = 4;
int marker_obstacle = 5;

// Current time (defined as global since needed in weak forms)
double TIME = 0;

// Boundary condition types for x-velocity
int xvel_bc_type(int marker) {
  if (marker == marker_right) return BC_NONE;
  else return BC_ESSENTIAL;
}

// Boundary condition values for x-velocity
scalar xvel_bc_value(int marker, double x, double y) {
  if (marker == marker_left) {
    // time-dependent inlet velocity (parabolic profile)
    double val_y = VEL_INLET * y*(H-y) / (H/2.)/(H/2.); //parabolic profile with peak VEL_INLET at y = H/2
    if (TIME <= STARTUP_TIME) return val_y * TIME/STARTUP_TIME;
    else return val_y;
  }
  else return 0;
}

// Boundary condition types for y-velocity
int yvel_bc_type(int marker) {
  if (marker == marker_right) return BC_NONE;
  else return BC_ESSENTIAL;
}

int p_bc_type(int marker)
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

  // load the mesh file
  Mesh basemesh, mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &basemesh); // master mesh
  basemesh.refine_all_elements();
  basemesh.refine_towards_boundary(5, 2, false);

  mesh.copy(&basemesh);

  // initialize shapesets and the cache
  H1ShapesetBeuchler h1_shapeset;
  PrecalcShapeset h1_pss(&h1_shapeset);
#ifdef PRESSURE_IN_L2
  L2Shapeset l2_shapeset;
  PrecalcShapeset l2_pss(&l2_shapeset);
#endif

  // spaces for velocities and pressure
  H1Space xvel_space(&mesh, &h1_shapeset);
  H1Space yvel_space(&mesh, &h1_shapeset);
#ifdef PRESSURE_IN_L2
  L2Space p_space(&mesh, &l2_shapeset);
#else
  H1Space p_space(&mesh, &h1_shapeset);
#endif

  // initialize boundary conditions
  xvel_space.set_bc_types(xvel_bc_type);
  xvel_space.set_bc_values(xvel_bc_value);
  yvel_space.set_bc_types(yvel_bc_type);
  p_space.set_bc_types(p_bc_type);

  // set velocity and pressure polynomial degrees
  xvel_space.set_uniform_order(P_INIT_VEL);
  yvel_space.set_uniform_order(P_INIT_VEL);
  p_space.set_uniform_order(P_INIT_PRESSURE);

  // assign degrees of freedom
  int ndof = assign_dofs(3, &xvel_space, &yvel_space, &p_space);

  // solutions for the Newton's iteration and time stepping
  Solution xvel_crs, yvel_crs, xvel_fine, yvel_fine, p_crs, p_fine;
  Solution xvel_prev_time, yvel_prev_time, xvel_prev_newton, yvel_prev_newton, p_prev;
  xvel_prev_time.set_zero(&mesh);
  yvel_prev_time.set_zero(&mesh);
  xvel_prev_newton.set_zero(&mesh);
  yvel_prev_newton.set_zero(&mesh);
  p_prev.set_zero(&mesh);

  // set up weak formulation
  WeakForm wf(3);
  if (NEWTON) {
    wf.add_biform(0, 0, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_biform(0, 0, callback(newton_bilinear_form_unsym_0_0), H2D_UNSYM, H2D_ANY, 2, &xvel_prev_newton, &yvel_prev_newton);
    wf.add_biform(0, 1, callback(newton_bilinear_form_unsym_0_1), H2D_UNSYM, H2D_ANY, 1, &xvel_prev_newton);
    wf.add_biform(0, 2, callback(bilinear_form_unsym_0_2), H2D_ANTISYM);
    wf.add_biform(1, 0, callback(newton_bilinear_form_unsym_1_0), H2D_UNSYM, H2D_ANY, 1, &yvel_prev_newton);
    wf.add_biform(1, 1, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_biform(1, 1, callback(newton_bilinear_form_unsym_1_1), H2D_UNSYM, H2D_ANY, 2, &xvel_prev_newton, &yvel_prev_newton);
    wf.add_biform(1, 2, callback(bilinear_form_unsym_1_2), H2D_ANTISYM);
    wf.add_liform(0, callback(newton_F_0), H2D_ANY, 5, &xvel_prev_time, &yvel_prev_time, &xvel_prev_newton, &yvel_prev_newton, &p_prev);
    wf.add_liform(1, callback(newton_F_1), H2D_ANY, 5, &xvel_prev_time, &yvel_prev_time, &xvel_prev_newton, &yvel_prev_newton, &p_prev);
    wf.add_liform(2, callback(newton_F_2), H2D_ANY, 2, &xvel_prev_newton, &yvel_prev_newton);
  }
  else {
    wf.add_biform(0, 0, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_biform(0, 0, callback(simple_bilinear_form_unsym_0_0_1_1), H2D_UNSYM, H2D_ANY, 2, &xvel_prev_time, &yvel_prev_time);
    wf.add_biform(1, 1, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_biform(1, 1, callback(simple_bilinear_form_unsym_0_0_1_1), H2D_UNSYM, H2D_ANY, 2, &xvel_prev_time, &yvel_prev_time);
    wf.add_biform(0, 2, callback(bilinear_form_unsym_0_2), H2D_ANTISYM);
    wf.add_biform(1, 2, callback(bilinear_form_unsym_1_2), H2D_ANTISYM);
    wf.add_liform(0, callback(simple_linear_form), H2D_ANY, 1, &xvel_prev_time);
    wf.add_liform(1, callback(simple_linear_form), H2D_ANY, 1, &yvel_prev_time);
  }

  // visualization
  VectorView vview("velocity [m/s]", 0, 0, 1500, 470);
  ScalarView pview("pressure [Pa]", 0, 530, 1500, 470);
  vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  //pview.set_min_max_range(-0.9, 1.0);
  pview.fix_scale_width(80);
  pview.show_mesh(true);

  // matrix solver
  UmfpackSolver umfpack;

  // linear system
  LinSystem ls(&wf, &umfpack);

  // nonlinear system
  NonlinSystem nls(&wf, &umfpack);

  if (NEWTON) {
    // set up the nonlinear system
    nls.set_spaces(3, &xvel_space, &yvel_space, &p_space);
#ifdef PRESSURE_IN_L2
    nls.set_pss(3, &h1_pss, &h1_pss, &l2_pss);
#else
    nls.set_pss(1, &h1_pss);
#endif
  }
  else {
    // set up the linear system
    ls.set_spaces(3, &xvel_space, &yvel_space, &p_space);
#ifdef PRESSURE_IN_L2
    ls.set_pss(3, &h1_pss, &h1_pss, &l2_pss);
#else
    ls.set_pss(1, &h1_pss);
#endif
  }

  // time-stepping loop
  char title[100];
  int num_time_steps = T_FINAL / TAU;
  for (int i = 1; i <= num_time_steps; i++)
  {
    TIME += TAU;
    info("---- Time step %d, time = %g:", i, TIME);

    // initial mesh and spaces
    mesh.copy(&basemesh);
    xvel_space.set_uniform_order(P_INIT_VEL);
    yvel_space.set_uniform_order(P_INIT_VEL);
    p_space.set_uniform_order(P_INIT_PRESSURE);

    // space adaptivity
    bool done = false; int at = 0;
    do
    {
      info("---- Adaptivity step %d:", at++);

      // assign degrees of freedom
      ndof = assign_dofs(3, &xvel_space, &yvel_space, &p_space);
      if (ndof >= NDOF_STOP) {
        done = true;
        break;
      }      

      if (NEWTON)
      {
        int it = 0; double res_l2_norm;
        do
        {
          info("---- Newton iter %d (coarse mesh):", it++);

          nls.assemble();
          nls.solve(3, &xvel_prev_newton, &yvel_prev_newton, &p_prev);
          res_l2_norm = nls.get_residuum_l2_norm();
          info("Residuum L2 norm: %g", res_l2_norm);
          if (it == 1) res_l2_norm = 100.0;
        }
        while (res_l2_norm > NEWTON_TOL && it <= NEWTON_MAX_ITER);

        sprintf(title, "Velocity, time %g s", TIME);
        vview.set_title(title);
        vview.show(&xvel_prev_newton, &yvel_prev_newton, H2D_EPS_LOW);
        sprintf(title, "Pressure, time %g s", TIME);
        pview.set_title(title);
        pview.show(&p_prev);

        xvel_crs.copy(&xvel_prev_newton);
        yvel_crs.copy(&yvel_prev_newton);

        it = 0;
        RefNonlinSystem refnls(&nls, 0);
        refnls.prepare();
        do
        {
          info("---- Newton iter %d (fine mesh):", it++);

          refnls.assemble();
          refnls.solve(3, &xvel_prev_newton, &yvel_prev_newton, &p_prev);
          res_l2_norm = refnls.get_residuum_l2_norm();
          info("Residuum L2 norm: %g", res_l2_norm);
        }
        while (res_l2_norm > NEWTON_TOL && it <= NEWTON_MAX_ITER);

        xvel_fine.copy(&xvel_prev_newton);
        yvel_fine.copy(&yvel_prev_newton);

      }
      else // linearized NS
      {
        ls.assemble();
        ls.solve(3, &xvel_crs, &yvel_crs, &p_crs);

        sprintf(title, "Velocity, time %g", TIME);
        vview.set_title(title);
        vview.show(&xvel_crs, &yvel_crs, H2D_EPS_LOW);
        sprintf(title, "Pressure, time %g", TIME);
        pview.set_title(title);
        pview.show(&p_crs);

        RefSystem refls(&ls, 0);
        refls.assemble();
        refls.solve(3, &xvel_fine, &yvel_fine, &p_fine);
      }

      DXDYFilter crs_mag(mag, &xvel_crs, &yvel_crs);
      DXDYFilter fine_mag(mag, &xvel_fine, &yvel_fine);
      double space_err = 100 * l2_error(&crs_mag, &fine_mag);
      info("Velocity rel error est %g%%", space_err);

      H1OrthoHP hp(1, &xvel_space);
      hp.set_biform(0, 0, callback(l2_form));
      space_err = hp.calc_error(&xvel_crs, &xvel_fine) * 100;
      info("L2 error (xvel) %g%%", space_err);
      if (space_err > SPACE_TOL) hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE, ISO_ONLY, MESH_REGULARITY, CONV_EXP);
      else done = true;

      xvel_space.set_uniform_order(P_INIT_VEL);
      yvel_space.set_uniform_order(P_INIT_VEL);
      p_space.set_uniform_order(P_INIT_PRESSURE);
    }
    while (!done);

    xvel_prev_time = xvel_fine;
    yvel_prev_time = yvel_fine;
  }

  // wait for keyboard or mouse input
  View::wait();
  return 0;
}
