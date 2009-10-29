#include "hermes2d.h"
#include "solver_umfpack.h"

#include "_hermes2d_api.h"

// The time-dependent laminar incompressible Navier-Stokes equations are
// discretized in time via the implicit Euler method. The convective term
// is linearized simply by replacing the velocity in front of the nabla
// operator with the velocity from last time step. Velocity is approximated
// using continuous elements, and pressure by means or discontinuous (L2)
// elements. This makes the velocity discretely divergence-free, meaning
// that the integral of div(v) over every element is zero. The problem
// has a steady symmetric solution which is unstable. Note that after some
// time (around t = 100), numerical errors induce oscillations. The
// approximation becomes unsteady and thus diverges from the exact solution.
// Interestingly, this happens even with a completely symmetric mesh.
//
// PDE: incompressible Navier-Stokes equations in the form
// \partial v / \partial t - \Delta v / Re + (v \cdot \nabla) v + \nabla p = 0,
// div v = 0
//
// BC: u_1 is a time-dependent constant and u_2 = 0 on Gamma_4 (inlet)
//     u_1 = u_2 = 0 on Gamma_1 (bottom), Gamma_3 (top) and Gamma_5 (obstacle)
//     "do nothing" on Gamma_2 (outlet)
//
// TODO: Implement Crank-Nicolson so that comparisons with implicit Euler can be made
//
// The following parameters can be played with:

const double RE = 1000.0;            // Reynolds number
const double VEL_INLET = 1.0;        // inlet velocity (reached after STARTUP_TIME)
const double STARTUP_TIME = 1.0;     // during this time, inlet velocity increases gradually
                                     // from 0 to VEL_INLET, then it stays constant
const double TAU = 0.5;              // time step
const double FINAL_TIME = 3000.0;    // length of time interval
const int P_INIT_VEL = 1;            // polynomial degree for velocity components
const int P_INIT_w0 = 1;       // polynomial degree for pressure
const int P_INIT_w4 = 1;       // polynomial degree for pressure
                                     // Note: P_INIT_VEL should always be greater than
                                     // P_INIT_PRESSURE because of the inf-sup condition
const double H = 5.0;                // domain height (necessary to define the parabolic
                                     // velocity profile at inlet)

const double R = 8.31;            // Gas constant
const double c_v = 1.0;            // specific heat capacity
const double g = 0;            // gravitational acceleration

//  boundary markers
int marker_bottom = 1;
int marker_right  = 2;
int marker_top = 3;
int marker_left = 4;
int marker_obstacle = 5;

// global time variable
double TIME = 0;

int bc_type(int marker) {
  return BC_ESSENTIAL;
}

scalar s0_bc_value(int marker, double x, double y) {
    return 1;
}

scalar s1_bc_value(int marker, double x, double y) {
    return 0.1;
}

scalar s3_bc_value(int marker, double x, double y) {
    return 0;
}

scalar s4_bc_value(int marker, double x, double y) {
    return 1;
}

scalar w0_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
    return 1;
}

scalar w1_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
    return 0.1;
}

scalar w3_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
    return 0;
}

scalar w4_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
    return 1;
}


double A_x(int i, int j, double w0, double w1, double w3, double w4)
{
    if (i == 0 && j == 0)
        return 0;
    else if (i == 0 && j == 1)
        return 1;
    else if (i == 0 && j == 2)
        return 0;
    else if (i == 0 && j == 3)
        return 0;

    else if (i == 1 && j == 0)
        return -w1*w1/(w0*w0) + R * (w1*w1 + w3*w3)/(2*c_v * w0*w0);
    else if (i == 1 && j == 1)
        return 2*w1/w0 - R * w1 / (c_v * w0);
    else if (i == 1 && j == 2)
        return -R * w3 / (c_v * w0);
    else if (i == 1 && j == 3)
        return R/c_v;

    else if (i == 2 && j == 0)
        return -w1*w3/(w0*w0);
    else if (i == 2 && j == 1)
        return w3/w0;
    else if (i == 2 && j == 2)
        return w1/w0;
    else if (i == 2 && j == 3)
        return 0;

    else if (i == 3 && j == 0)
        return -w1*w4/(w0*w0) - w1/(w0*w0) * R/c_v
            * (w4 - (w1*w1+w3*w3)/(2*w0)) + w1/w0 * R/c_v
            * (w1*w1+w3*w3)/(2*w0*w0);
        // or equivalently:
        //return w1/w0 * (R/c_v * (w1*w1+w3*w3)/(w0*w0) - (R/c_v + 1) * w4/w0);
    else if (i == 3 && j == 1)
        return w4/w0 + 1/w0 * R/c_v
            * (w4 - (w1*w1+w3*w3)/(2*w0)) - R/c_v
            * w1*w1/(w0*w0);
    else if (i == 3 && j == 2)
        return - R/c_v * w1*w3/(w0*w0);
    else if (i == 3 && j == 3)
        return w1/w0 + R/c_v * w1/w0;

    printf("i=%d, j=%d;\n", i, j);
    error("Invalid index.");
}

double A_z(int i, int j, double w0, double w1, double w3, double w4)
{
    if (i == 0 && j == 0)
        return 0;
    else if (i == 0 && j == 1)
        return 0;
    else if (i == 0 && j == 2)
        return 1;
    else if (i == 0 && j == 3)
        return 0;

    else if (i == 1 && j == 0)
        return -w3*w1/(w0*w0);
    else if (i == 1 && j == 1)
        return w3/w0;
    else if (i == 1 && j == 2)
        return w1/w0;
    else if (i == 1 && j == 3)
        return 0;

    else if (i == 2 && j == 0)
        return -w3*w3/(w0*w0) + R * (w1*w1 + w3*w3)/(2*c_v * w0*w0);
    else if (i == 2 && j == 1)
        return -R * w1 / (c_v * w0);
    else if (i == 2 && j == 2)
        return 2*w3/w0 - R * w3 / (c_v * w0);
    else if (i == 2 && j == 3)
        return R/c_v;

    else if (i == 3 && j == 0)
        return -w3*w4/(w0*w0) - w3/(w0*w0) * R/c_v
            * (w4 - (w1*w1+w3*w3)/(2*w0)) + w3/w0 * R/c_v
            * (w1*w1+w3*w3)/(2*w0*w0);
        // or equivalently:
        //return w1/w0 * (R/c_v * (w1*w1+w3*w3)/(w0*w0) - (R/c_v + 1) * w4/w0);
    else if (i == 3 && j == 1)
        return - R/c_v * w3*w1/(w0*w0);
    else if (i == 3 && j == 2)
        return w4/w0 + 1/w0 * R/c_v
            * (w4 - (w1*w1+w3*w3)/(2*w0)) - R/c_v
            * w3*w3/(w0*w0);
    else if (i == 3 && j == 3)
        return w3/w0 + R/c_v * w3/w0;

    error("Invalid index.");
}

double f_x(int i, double w0, double w1, double w3, double w4)
{
    if (i == 0)
        return w1;
    else if (i == 1)
        return w1*w1/w0 + R/c_v * (w4 - (w1*w1+w3*w3)/(2*w0));
    else if (i == 2)
        return w1*w3/w0;
    else if (i == 3)
        return w1/w0 * (w4 + R/c_v * (w4 - (w1*w1+w3*w3)/(2*w0)));

    error("Invalid index.");
}

double f_z(int i, double w0, double w1, double w3, double w4)
{
    if (i == 0)
        return w3;
    else if (i == 1)
        return w3*w1/w0;
    else if (i == 2)
        return w3*w3/w0 + R/c_v * (w4 - (w1*w1+w3*w3)/(2*w0));
    else if (i == 3)
        return w3/w0 * (w4 + R/c_v * (w4 - (w1*w1+w3*w3)/(2*w0)));

    error("Invalid index.");
}

double w(int i, double w0, double w1, double w3, double w4)
{
    if (i == 0)
        return w0;
    else if (i == 1)
        return w1;
    else if (i == 2)
        return w3;
    else if (i == 3)
        return w4;

    error("Invalid index.");
}

double test1(double w0, double w1, double w3, double w4)
{
    for (int i=0; i<4; i++) {
        double r=0;
        for (int j=0; j<4; j++) {
            r += A_x(i, j, w0, w1, w3, w4) * w(j, w0, w1, w3, w4);
        }
        r -= f_x(i, w0, w1, w3, w4);
        printf("result: %d: %f\n", i, r);
    }
}

double test3(double w0, double w1, double w3, double w4)
{
    for (int i=0; i<4; i++) {
        double r=0;
        for (int j=0; j<4; j++) {
            r += A_z(i, j, w0, w1, w3, w4) * w(j, w0, w1, w3, w4);
        }
        r -= f_z(i, w0, w1, w3, w4);
        printf("result: %d: %f\n", i, r);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar B_ij(int _i, int _j, int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    double delta_ij;
    if (_i == _j)
        delta_ij = 1;
    else
        delta_ij = 0;
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (
                (u->val[i] * v->val[i]) * delta_ij / TAU +
                -A_x(_i, _j, ext->fn[0]->val[i],
                    ext->fn[1]->val[i],
                    ext->fn[2]->val[i],
                    ext->fn[3]->val[i])
                * u->val[i] * v->dx[i] +
                -A_z(_i, _j, ext->fn[0]->val[i],
                    ext->fn[1]->val[i],
                    ext->fn[2]->val[i],
                    ext->fn[3]->val[i])
                * u->val[i] * v->dy[i]
                );
    return result;
}

template<typename Real, typename Scalar>
Scalar S_ij(int _i, int _j, int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (
                A_x(_i, _j, ext->fn[0]->val[i],
                    ext->fn[1]->val[i],
                    ext->fn[2]->val[i],
                    ext->fn[3]->val[i])
                * u->val[i] * v->val[i] +
                A_z(_i, _j, ext->fn[0]->val[i],
                    ext->fn[1]->val[i],
                    ext->fn[2]->val[i],
                    ext->fn[3]->val[i])
                * u->val[i] * v->val[i]
                );
    return result;
}

template<typename Real, typename Scalar>
Scalar B_00(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(0, 0, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_01(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(0, 1, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_02(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(0, 2, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_03(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(0, 3, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_10(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(1, 0, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_11(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(1, 1, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_12(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(1, 2, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_13(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(1, 3, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_20(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(2, 0, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_21(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(2, 1, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_22(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(2, 2, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_23(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(2, 3, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_30(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(3, 0, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_31(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(3, 1, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_32(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(3, 2, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar B_33(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return B_ij(3, 3, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar l_0(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (ext->fn[0]->val[i] * v->val[i]) / TAU;
    return result;
}

template<typename Real, typename Scalar>
Scalar l_1(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (ext->fn[1]->val[i] * v->val[i]) / TAU;
    return result;
}

template<typename Real, typename Scalar>
Scalar l_2(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (ext->fn[2]->val[i]/TAU + ext->fn[0]->val[i]*g) * \
                  v->val[i] / TAU;
    return result;
}

template<typename Real, typename Scalar>
Scalar l_3(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (ext->fn[3]->val[i] * v->val[i]) / TAU;
    return result;
}

template<typename Real, typename Scalar>
Scalar S_00(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(0, 0, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_01(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(0, 1, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_02(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(0, 2, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_03(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(0, 3, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_10(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(1, 0, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_11(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(1, 1, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_12(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(1, 2, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_13(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(1, 3, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_20(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(2, 0, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_21(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(2, 1, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_22(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(2, 2, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_23(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(2, 3, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_30(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(3, 0, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_31(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(3, 1, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_32(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(3, 2, n, wt, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar S_33(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return S_ij(3, 3, n, wt, u, v, e, ext);
}

Ord B_order(int n, double *wt, Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
     return Ord(20);
}

#define callback_bf(a) a<double, scalar>, B_order

int main(int argc, char* argv[])
{
  // import hermes2d

  printf("Importing hermes1d\n");
  // Initialize Python
  Py_Initialize();
  PySys_SetArgv(argc, argv);
  if (import_hermes2d___hermes2d())
      throw std::runtime_error("hermes2d failed to import.");
  cmd("print 'Python initialized'");

  // load the mesh file
  Mesh mesh;
  mesh.load("domain-quad.mesh"); // unstructured triangular mesh available in domain-tri.mesh

  // a-priori mesh refinements
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  //mesh.refine_towards_boundary(5, 4, false);
  //mesh.refine_towards_boundary(1, 4);
  //mesh.refine_towards_boundary(3, 4);

  // display the mesh
  //MeshView mview("Navier-Stokes Example - Mesh", 100, 100, 1100, 400);
  //mview.show(&mesh);
  //mview.wait_for_keypress();

  // initialize the shapesets and the cache
  H1Shapeset shapeset_h1;
  PrecalcShapeset pss_h1(&shapeset_h1);

  // H1 spaces for velocities and L2 for pressure
  H1Space s0(&mesh, &shapeset_h1);
  H1Space s1(&mesh, &shapeset_h1);
  H1Space s3(&mesh, &shapeset_h1);
  H1Space s4(&mesh, &shapeset_h1);

  // initialize boundary conditions
  s0.set_bc_types(bc_type);
  s0.set_bc_values(s0_bc_value);
  s1.set_bc_types(bc_type);
  s1.set_bc_values(s1_bc_value);
  s3.set_bc_types(bc_type);
  s3.set_bc_values(s3_bc_value);
  s4.set_bc_types(bc_type);
  s4.set_bc_values(s4_bc_value);

  // set velocity and pressure polynomial degrees
  s0.set_uniform_order(P_INIT_w0);
  s1.set_uniform_order(P_INIT_VEL);
  s3.set_uniform_order(P_INIT_VEL);
  s4.set_uniform_order(P_INIT_w4);

  // assign degrees of freedom
  int ndofs = 0;
  ndofs += s0.assign_dofs(ndofs);
  ndofs += s1.assign_dofs(ndofs);
  ndofs += s3.assign_dofs(ndofs);
  ndofs += s4.assign_dofs(ndofs);

  // initial condition: xprev and yprev are zero
  Solution w0_prev, w1_prev, w3_prev, w4_prev;
  w0_prev.set_exact(&mesh, w0_init);
  w1_prev.set_exact(&mesh, w1_init);
  w3_prev.set_exact(&mesh, w3_init);
  w4_prev.set_exact(&mesh, w4_init);

  // set up weak formulation
  WeakForm wf(4);
  wf.add_biform(0, 0, callback_bf(B_00), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(0, 1, callback_bf(B_01), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(0, 2, callback_bf(B_02), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(0, 3, callback_bf(B_03), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(1, 0, callback_bf(B_10), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(1, 1, callback_bf(B_11), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(1, 2, callback_bf(B_12), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(1, 3, callback_bf(B_13), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(2, 0, callback_bf(B_20), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(2, 1, callback_bf(B_21), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(2, 2, callback_bf(B_22), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(2, 3, callback_bf(B_23), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(3, 0, callback_bf(B_30), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(3, 1, callback_bf(B_31), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(3, 2, callback_bf(B_32), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_biform(3, 3, callback_bf(B_33), UNSYM, ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);

  wf.add_liform(0, callback(l_0), ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_liform(1, callback(l_1), ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_liform(2, callback(l_2), ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);
  wf.add_liform(3, callback(l_3), ANY, 4, &w0_prev, &w1_prev,
          &w3_prev, &w4_prev);


  // visualization
  VectorView w13_view("Current Density [m/s]", 0, 0, 1500, 470);
  ScalarView w0_view("Mass Density [Pa]", 0, 530, 1500, 470);
  ScalarView w4_view("Energy [Pa]", 0, 530, 1500, 470);
  //w13_view.set_min_max_range(0, 1.6);
  w0_view.show_mesh(false);
  w4_view.show_mesh(false);
  // fixing scale width (for nicer videos). Note: creation of videos is
  // discussed in a separate example
  //vview.fix_scale_width(5);
  //pview.fix_scale_width(5);

  // set up the linear system
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(4, &s0, &s1, &s3, &s4);
  sys.set_pss(4, &pss_h1, &pss_h1, &pss_h1, &pss_h1);
  //sys.set_pss(1, &pss_h1);

  // main loop
  char title[100];
  int num_time_steps = FINAL_TIME / TAU;
  for (int i = 1; i <= num_time_steps; i++)
  {
    TIME += TAU;

    cmd("print 'Iteration'");
    info("\n---- Time step %d, time = %g -----------------------------------", i, TIME);

    // this is needed to update the time-dependent boundary conditions
    ndofs = 0;
    ndofs += s0.assign_dofs(ndofs);
    ndofs += s1.assign_dofs(ndofs);
    ndofs += s3.assign_dofs(ndofs);
    ndofs += s4.assign_dofs(ndofs);

    // assemble and solve
    Solution w0_sln, w1_sln, w3_sln, w4_sln;
    sys.assemble();
    sys.solve(4, &w0_sln, &w1_sln, &w3_sln, &w4_sln);

    // visualization
    sprintf(title, "Current density, time %g", TIME);
    w13_view.set_title(title);
    w13_view.show(&w1_prev, &w3_prev, EPS_LOW);
    sprintf(title, "Mass Density, time %g", TIME);
    w0_view.set_title(title);
    w0_view.show(&w0_sln);
    sprintf(title, "Energy, time %g", TIME);
    w4_view.set_title(title);
    w4_view.show(&w0_sln);

    w0_prev = w0_sln;
    w1_prev = w1_sln;
    w3_prev = w3_sln;
    w4_prev = w4_sln;
  }

  View::wait();
}
