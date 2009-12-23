#include "hermes2d.h"
#include "solver_umfpack.h"

// This test makes sure that example 04-bc-dirichlet works correctly.
// CAUTION: This test will fail when any changes to the shapeset
// are made, but it is easy to fix (see below).

double CONST_F = -4.0;       // constant right-hand side

// boundary condition type (essential = Dirichlet)
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// function values for Dirichlet boundary markers
scalar bc_values(int marker, double x, double y)
{
  return (-CONST_F/4.0)*(x*x + y*y);
}

// return the value \int \nabla u . \nabla v dx
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

// return the value \int v dx
template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return CONST_F*int_v<Real, Scalar>(n, wt, v);
}

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");
  mesh.refine_all_elements();

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form));
  wf.add_liform(0, callback(linear_form));

  // initialize the linear system and solver
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);

  // testing n_dof and correctness of solution vector 
  // for p_init = 1, 2, ..., 10
  int success = 1;
  for (int p_init = 1; p_init <= 10; p_init++) {
    printf("********* p_init = %d *********\n", p_init);
    space.set_uniform_order(p_init);
    space.assign_dofs();

    // assemble the stiffness matrix and solve the system
    Solution sln;
    sys.assemble();
    sys.solve(1, &sln);

    scalar *sol_vector;
    int n_dof;
    sys.get_solution_vector(sol_vector, n_dof);
    printf("n_dof = %d\n", n_dof);
    double sum = 0;
    for (int i=0; i < n_dof; i++) sum += sol_vector[i];
    printf("coefficient sum = %g\n", sum);

    // Actual test. The values of 'sum' depend on the 
    // current shapeset. If you change the shapeset, 
    // you need to correct these numbers. 
    if (p_init == 1 && fabs(sum - 1.7251) > 1e-3) success = 0;
    if (p_init == 2 && fabs(sum - 3.79195) > 1e-3) success = 0;
    if (p_init == 3 && fabs(sum - 3.80206) > 1e-3) success = 0;
    if (p_init == 4 && fabs(sum - 3.80156) > 1e-3) success = 0;
    if (p_init == 5 && fabs(sum - 3.80155) > 1e-3) success = 0;
    if (p_init == 6 && fabs(sum - 3.80154) > 1e-3) success = 0;
    if (p_init == 7 && fabs(sum - 3.80154) > 1e-3) success = 0;
    if (p_init == 8 && fabs(sum - 3.80153) > 1e-3) success = 0;
    if (p_init == 9 && fabs(sum - 3.80152) > 1e-3) success = 0;
    if (p_init == 10 && fabs(sum - 3.80152) > 1e-3) success = 0;
  }

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
