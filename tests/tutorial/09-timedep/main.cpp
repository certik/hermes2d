#include "hermes2d.h"

// This test makes sure that example 09-timedep works correctly.
// CAUTION: This test will fail when any changes to the shapeset
// are made, but it is easy to fix (see below).

const int P_INIT = 3;            // initial polynomial degree in elements
const int INIT_REF_NUM = 1;      // number of initial uniform refinements
const double TAU = 200.0;        // time step in seconds

// Problem constants
const double T_INIT = 10;        // temperature of the ground (also initial temperature)
const double ALPHA = 10;         // heat flux coefficient for Newton's boundary condition
const double LAMBDA = 1e5;       // thermal conductivity of the material
const double HEATCAP = 1e6;      // heat capacity
const double RHO = 3000;         // material density
const double FINAL_TIME = 2100;  // length of time interval (24 hours) in seconds

// Global variable.
double TIME = 0;

// Time-dependent exterior temperature.
double temp_ext(double t) {
  return T_INIT + 10. * sin(2*M_PI*t/FINAL_TIME);
}

// Boundary markers.
int marker_ground = 1;
int marker_air = 2;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == marker_ground) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Function values for Dirichlet boundary markers.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return T_INIT;
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return HEATCAP * RHO * int_u_v<Real, Scalar>(n, wt, u, v) / TAU +
         LAMBDA * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return HEATCAP * RHO * int_u_v<Real, Scalar>(n, wt, ext->fn[0], v) / TAU;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return LAMBDA * ALPHA * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return LAMBDA * ALPHA * temp_ext(TIME) * int_v<Real, Scalar>(n, wt, v);
}


int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("cathedral.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(2, 5);

  // Initialize an H1 space.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Set initial condition.
  Solution tsln;
  tsln.set_const(&mesh, T_INIT);

  // Initialize weak formulation.
  WeakForm wf;
  wf.add_matrix_form(bilinear_form<double, double>, bilinear_form<Ord, Ord>);
  wf.add_matrix_form_surf(bilinear_form_surf<double, double>, bilinear_form_surf<Ord, Ord>, marker_air);
  wf.add_vector_form(linear_form<double, double>, linear_form<Ord, Ord>, H2D_ANY, &tsln);
  wf.add_vector_form_surf(linear_form_surf<double, double>, linear_form_surf<Ord, Ord>, marker_air);

  // Initialize linear system.
  LinSystem ls(&wf, &space);

  // time stepping
  int nsteps = (int)(FINAL_TIME/TAU + 0.5);
  bool rhsonly = false;
  for(int n = 1; n <= nsteps; n++)
  {
    info("---- Time step %d, time %3.5f, ext_temp %g", ts, TIME, temp_ext(TIME));

    // Assemble and solve.
    ls.assemble(rhsonly);
    rhsonly = true;
    ls.solve(&tsln);

    // Update the time variable.
    TIME += TAU;
  }

  scalar *sol_vector;
  int n_dof;
  ls.get_solution_vector(sol_vector, n_dof);
  printf("n_dof = %d\n", n_dof);
  double sum = 0;
  for (int i=0; i < n_dof; i++) sum += sol_vector[i];
  printf("coefficient sum = %g\n", sum);

  // Actual test. The value of 'sum' depend on the
  // current shapeset. If you change the shapeset,
  // you need to correct this number.
  int success = 1;
  if (fabs(sum - 9122.66) > 1e-1) success = 0;

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  if (success == 1) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}

