#include "forms.h"

#include "_hermes2d_api.h"

const double TAU = 1;  // this is in seconds

const double R = 287.14;            // Gas constant
//const double T2 = 310;
//const double T1 = 275;
//const double T1 = 310;
const double T_0 = 300.5;
const double p_0 = 100000;
const double rho_0 = p_0/(R*T_0);
const double c_v = 20.8;            // specific heat capacity
const double g = 9.81;            // gravitational acceleration (set to 0 for now)

//  boundary markers
#define marker_bottom 1
#define marker_right 2
#define marker_top 3
#define marker_left 4

int s0_bc_type(int marker) {
    return BC_NATURAL;
}

int s1_bc_type(int marker) {
    if (marker == marker_top || marker == marker_bottom)
        return BC_NATURAL;
    else
        return BC_ESSENTIAL;
}

int s3_bc_type(int marker) {
    if (marker == marker_left || marker == marker_right)
        return BC_NATURAL;
    else
        return BC_ESSENTIAL;
}

int s4_bc_type(int marker) {
    if (marker=marker_bottom)
        return BC_ESSENTIAL;
    else
        return BC_NATURAL;
}

/*
scalar s0_bc_value(int marker, double x, double y) {
    return 1;
}
*/

scalar s1_bc_value(int marker, double x, double y) {
    return 0;
}

scalar s3_bc_value(int marker, double x, double y) {
    return 0;
}

scalar s4_bc_value(int marker, double x, double y) {
    return rho_0 * c_v * T_0 * (1 + tanh(x/50000));
}

// Empirical relations of initial distributions valid for 0 <= z <= 10km
// "z" in p_z and T_z is in "km", so don't forget to convert it from meters
// temperature T in Kelvin
#define T_z(z) (T_0 - 8.3194*(z) + 0.2932*(z)*(z) - 0.0109*(z)*(z)*(z))
// presure p in Pascals
#define p_z(z) (p_0 - 11476*(z) + 529.54*(z)*(z) - 9.38*(z)*(z)*(z))

#define theta_z(z) (T_z(z) * pow(p_0/p_z(z), 0.287))
#define gamma 1.4
#define rho_z(z) (p_0 / (R * theta_z(z)) * pow(p_z(z)/p_0, 1./gamma))
#define rho_init(x, y) (rho_z(y/1000.))
#define T_init(x, y) (T_z(y/1000.))

scalar w0_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
    return rho_init(x, y);
}

scalar w1_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
    return 0;
}

scalar w3_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
    return 0;
}

scalar w4_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
//    printf("y=%f, rho=%f\n", y, rho_init(x, y));
    return rho_init(x, y) * T_init(x, y) * c_v;
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
        result += wt[i] * (ext->fn[2]->val[i]/TAU - ext->fn[0]->val[i]*g) * \
                  v->val[i];
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

void register_bc(H1Space &s0, H1Space &s1, H1Space &s3, H1Space &s4)
{
    s0.set_bc_types(s0_bc_type);
    //s0.set_bc_values(s0_bc_value);
    s1.set_bc_types(s1_bc_type);
    s1.set_bc_values(s1_bc_value);
    s3.set_bc_types(s3_bc_type);
    s3.set_bc_values(s3_bc_value);
    s4.set_bc_types(s4_bc_type);
    s4.set_bc_values(s4_bc_value);
}

void set_ic(Mesh &mesh, Solution &w0, Solution &w1, Solution &w3, Solution &w4)
{
    w0.set_exact(&mesh, w0_init);
    w1.set_exact(&mesh, w1_init);
    w3.set_exact(&mesh, w3_init);
    w4.set_exact(&mesh, w4_init);
}

void register_forms(WeakForm &wf, Solution &w0_prev, Solution &w1_prev,
        Solution &w3_prev, Solution &w4_prev)
{
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

    wf.add_liform(0, callback(l_0), ANY, 4, &w0_prev, &w1_prev, &w3_prev,
            &w4_prev);
    wf.add_liform(1, callback(l_1), ANY, 4, &w0_prev, &w1_prev, &w3_prev,
            &w4_prev);
    wf.add_liform(2, callback(l_2), ANY, 4, &w0_prev, &w1_prev, &w3_prev,
            &w4_prev);
    wf.add_liform(3, callback(l_3), ANY, 4, &w0_prev, &w1_prev, &w3_prev,
            &w4_prev);
}
