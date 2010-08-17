#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function.h"

// This test makes sure that example 18-newton-timedep-heat works correctly.

const int INIT_GLOB_REF_NUM = 3;       // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 4;        // Number of initial refinements towards boundary.
const int P_INIT = 2;                  // Initial polynomial degree.
const double TAU = 0.2;                // Time step.
const double T_FINAL = 5.0;            // Time interval length.
const double NEWTON_TOL = 1e-6;        // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;       // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Thermal conductivity (temperature-dependent).
// Note: for any u, this function has to be positive.
template<typename Real>
Real lam(Real u)
{
  return 1 + pow(u, 4);
}

// Derivative of the thermal conductivity with respect to 'u'.
template<typename Real>
Real dlam_du(Real u) {
  return 4*pow(u, 3);
}

// This function is used to define Dirichlet boundary conditions.
double dir_lift(double x, double y, double& dx, double& dy) {
  dx = (y+10)/10.;
  dy = (x+10)/10.;
  return (x+10)*(y+10)/100.;
}

// Initial condition. It will be projected on the FE mesh 
// to obtain initial coefficient vector for the Newton's method.
scalar init_cond(double x, double y, double& dx, double& dy)
{
  return dir_lift(x, y, dx, dy);
}

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition markers.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  double dx, dy;
  return dir_lift(x, y, dx, dy);
}

// Heat sources (can be a general function of 'x' and 'y').
template<typename Real>
Real heat_src(Real x, Real y)
{
  return 1.0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);

  // Create an H1 space with default shapeset.
  H1Space* space = new H1Space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(space);
  info("ndof = %d.", ndof);

  // Solutions for the time stepping and the Newton's method.
  Solution u_prev_time;

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(jac), H2D_UNSYM, H2D_ANY);
  wf.add_vector_form(callback(res), H2D_ANY, &u_prev_time);

  // Project the initial condition on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  Vector* coeff_vec = new AVector(ndof); 
  info("Projecting initial condition to obtain initial vector for the Newton'w method.");
  project_global(space, H2D_H1_NORM, init_cond, &u_prev_time, coeff_vec);

  // Time stepping loop:
  double current_time = 0.0;
  int ts = 1;
  do {
    info("---- Time step %d, t = %g s.", ts, current_time); ts++;

    // Newton's method.
    info("Performing Newton's method.");
    bool verbose = true; // Default is false.
    if (!solve_newton(space, &wf, coeff_vec, matrix_solver, NEWTON_TOL, NEWTON_MAX_ITER, verbose))
      error("Newton's method did not converge.");

    // Update previous time level solution.
    u_prev_time.set_fe_solution(space, coeff_vec);

    // Update time.
    current_time += TAU;

  } while (current_time < T_FINAL);

  ndof = get_num_dofs(space);

#define ERROR_SUCCESS                                0
#define ERROR_FAILURE                               -1
  printf("ndof allowed = %d\n", 1000);
  printf("ndof actual = %d\n", ndof);
  if (ndof < 1000) {      // ndofs was 961 at the time this test was created.
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}

