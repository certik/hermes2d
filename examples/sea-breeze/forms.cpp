#include "forms.h"
#include "numerical_flux.h"

#include "_hermes2d_api.h"

const double TAU = 0.01/t_r;  // this is in seconds

const double T_0 = T_z(0);
const double p_0 = p_z(0);
const double rho_0 = rho_z(0);

int s0_bc_type(int marker) {
    return BC_NATURAL;
}

int s1_bc_type(int marker) {
    return BC_NATURAL;
}

int s3_bc_type(int marker) {
    return BC_NATURAL;
}

int s4_bc_type(int marker) {
    return BC_NATURAL;
}

// XXX: This is a hack, it should be made iteration independent
int iterations = 0;
void set_iteration(int i) {
    iterations = i;
}

#define rho_init(x, y) (rho_z(y*l_r))
#define T_init(x, y) (T_z(y*l_r))

double w0_init_num;
double w1_init_num;
double w3_init_num;
double w4_init_num;
double p_init_num;


scalar w0_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
    //w0_init_num = rho_z(0)/rho_r;
    w0_init_num = 1;
    return w0_init_num;
}

scalar w1_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
    //w1_init_num = rho_z(0)/rho_r * (20/u_r);
    w1_init_num = 1;
    return w1_init_num;
}

scalar w3_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
    //w3_init_num = rho_z(0)/rho_r * (0/u_r);
    w3_init_num = 0;
    return w3_init_num;
}


scalar w4_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
    //w4_init_num = rho_z(0) * T_z(0) * c_v / E_r;
    w4_init_num = 1;
    p_init_num = R/c_v * (w4_init_num - (w1_init_num*w1_init_num +
                w3_init_num*w3_init_num)/(2*w0_init_num));
    return w4_init_num;
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
    assert(u->val != NULL);
    assert(v->dx != NULL);
    assert(v->dy != NULL);
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
    double w0, w1, w3, w4;

    if (e->marker == marker_left) {
        printf("BC left: (%d, %d; x=%f y=%f)\n", _i, _j, e->x[0], e->y[0]);
        //printf("    state: (%f, %f, %f, %f)\n", w0, w1, w3, w4);
            printf("  vvv: %d\n", n);
            printf("   u:");
            for (int j = 0; j<n;j++)
                printf("%f ", u->val[j]);
            printf("\n");
            printf("   v:");
            for (int j = 0; j<n;j++)
                printf("%f ", v->val[j]);
            printf("\n");
    }
    Scalar result = 0;
    for (int i = 0; i < n; i++) {
        w0 = ext->fn[0]->val[i];
        w1 = ext->fn[1]->val[i];
        w3 = ext->fn[2]->val[i];
        w4 = ext->fn[3]->val[i];
        result += wt[i] * (
                A_x(_i, _j, w0, w1, w3, w4) * e->nx[i]
                +
                A_z(_i, _j, w0, w1, w3, w4) * e->ny[i]
                ) * u->val[i] * v->val[i];
    }
    if (e->marker == marker_left) {
        printf("   result: %f\n", result);
    }
    return result;
}

template<typename Real, typename Scalar>
Scalar s_i(int _i, int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    double w0, w1, w3, w4;
    Scalar result = 0;
    //printf("BC: n=%d; marker=%d\n", n, e->marker);
    for (int i = 0; i < n; i++) {
        w0 = ext->fn[0]->val[i];
        w1 = ext->fn[1]->val[i];
        w3 = ext->fn[2]->val[i];
        w4 = ext->fn[3]->val[i];
        double _rho = w0;
        double _u = w1/w0;
        double _w = w3/w0;
        double _E = w4;
        double _v2 = _u*_u+_w*_w;
        double _p = (kappa-1)*(_E - _rho*_v2/2);
        double _c = sqrt(kappa*_p/_rho);
        double _M = sqrt(_v2)/_c;
        double un = _u*e->nx[i] + _w*e->ny[i];
        // we build a new state w?_new that we want to impose
        double w0_new, w1_new, w3_new, w4_new;
        int right_or_left = 0;
        if (e->marker == marker_left || e->marker == marker_right)
            right_or_left = 1;
        else
            right_or_left = 0;
        right_or_left = 1;

        if  (!right_or_left) {
            _u = _u-un * e->nx[i];
            _w = _w-un * e->ny[i];
            w0_new = w0;
            w1_new = _u * w0;
            w3_new = _w * w0;
            w4_new = w4;
        } else {
            if (un > 0) {
                // outlet
                if (_M >= 1.) {
                    // supersonic
                    // we take everything from inside
                    w0_new = w0;
                    w1_new = w1;
                    w3_new = w3;
                    w4_new = w4;
                }
                else {
                    // subsonic
                    // take p from the outside state, the rest from the inside
                    w0_new = w0;
                    w1_new = w1;
                    w3_new = w3;
                    w4_new = p_init_num * c_v / R + (w1*w1+w3*w3)/(2*w0);
                }
            }
            else {
                // inlet
                if (_M >= 1.) {
                    // supersonic
                    // take everything from outside
                    w0_new = w0_init_num;
                    w1_new = w1_init_num;
                    w3_new = w3_init_num;
                    w4_new = w4_init_num;
                }
                else {
                    // subsonic
                    // take rho, u, w from the outside state
                    w0_new = w0_init_num;
                    w1_new = w1_init_num;
                    w3_new = w3_init_num;
                    // calculate E:
                    double v2 = w1_new*w1_new+w3_new*w3_new;
                    w4_new = _p * c_v / R + v2 / (2*w0_new);
                }
            }
        }
        // we only impose the difference to the inside state:
        w0_new -= w0;
        w1_new -= w1;
        w3_new -= w3;
        w4_new -= w4;
        result += wt[i] * (
                A_x(_i, 0, w0, w1, w3, w4) * w0_new * e->nx[i]
                +
                A_x(_i, 1, w0, w1, w3, w4) * w1_new * e->nx[i]
                +
                A_x(_i, 2, w0, w1, w3, w4) * w3_new * e->nx[i]
                +
                A_x(_i, 3, w0, w1, w3, w4) * w4_new * e->nx[i]
                +
                A_z(_i, 0, w0, w1, w3, w4) * w0_new * e->ny[i]
                +
                A_z(_i, 1, w0, w1, w3, w4) * w1_new * e->ny[i]
                +
                A_z(_i, 2, w0, w1, w3, w4) * w3_new * e->ny[i]
                +
                A_z(_i, 3, w0, w1, w3, w4) * w4_new * e->ny[i]
                ) * v->val[i];
    }
    if (e->marker == marker_left) {
        printf("BC left: (%d; x=%f y=%f) %f\n", _i, e->x[0], e->y[0],
                result);
        printf("    state: (%f, %f, %f, %f)\n", w0, w1, w3, w4);
    }
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
        result += wt[i] * (ext->fn[2]->val[i]/TAU - ext->fn[0]->val[i]*g/g_r) *
                  v->val[i];
    return result;
}

template<typename Real, typename Scalar>
Scalar l_3(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    insert_object("fn", array_double_c2numpy(ext->fn[3]->val, n));
    insert_object("v", array_double_c2numpy(v->val, n));
    insert_object("tau", double_c2py(TAU));
    cmd("print fn");
    cmd("print v");
    cmd("print tau");
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (ext->fn[3]->val[i] * v->val[i]) / TAU;
    return result;
}

template<typename Real, typename Scalar>
Scalar l_ord(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (ext->fn[3]->val[i] * v->val[i]) / TAU;
    return result;
}

template<typename Real, typename Scalar>
Scalar s_1(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return s_i(1, n, wt, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar s_2(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return s_i(2, n, wt, v, e, ext);
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

Ord L_order(int n, double *wt, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
     return Ord(20);
}

void register_bc(Space &s0, Space &s1, Space &s3, Space &s4)
{
    s0.set_bc_types(s0_bc_type);
    s1.set_bc_types(s1_bc_type);
    s3.set_bc_types(s3_bc_type);
    s4.set_bc_types(s4_bc_type);
}

void set_ic(Mesh &mesh, Solution &w0, Solution &w1, Solution &w3, Solution &w4)
{
    w0.set_exact(&mesh, w0_init);
    w1.set_exact(&mesh, w1_init);
    w3.set_exact(&mesh, w3_init);
    w4.set_exact(&mesh, w4_init);
}

#define callback_bf(a) a<double, scalar>, B_order

#define callback_lf_s(a) a<double, scalar>, L_order

#define callback_lf(a) a<double, scalar>, l_ord<Ord, Ord>

#define ADD_BF(i, j) wf.add_biform(i, j, callback_bf(B_##i##j), UNSYM, ANY, 4, &w0_prev, &w1_prev, &w3_prev, &w4_prev)

#define ADD_BF_S(i, j) wf.add_biform_surf(i, j, callback_bf(S_##i##j), ANY, 4, &w0_prev, &w1_prev, &w3_prev, &w4_prev)

#define ADD_LF(i) wf.add_liform(i, callback_lf(l_##i), ANY, 4, &w0_prev, &w1_prev, &w3_prev, &w4_prev);

#define ADD_LF_S(i) wf.add_liform_surf(i, callback_lf_s(s_##i), ANY, 4, &w0_prev, &w1_prev, &w3_prev, &w4_prev);

void register_forms(WeakForm &wf, Solution &w0_prev, Solution &w1_prev,
        Solution &w3_prev, Solution &w4_prev)
{
    ADD_BF(0, 0); ADD_BF(0, 1); ADD_BF(0, 2); ADD_BF(0, 3);
    ADD_BF(1, 0); ADD_BF(1, 1); ADD_BF(1, 2); ADD_BF(1, 3);
    ADD_BF(2, 0); ADD_BF(2, 1); ADD_BF(2, 2); ADD_BF(2, 3);
    ADD_BF(3, 0); ADD_BF(3, 1); ADD_BF(3, 2); ADD_BF(3, 3);

    ADD_BF_S(0, 0); ADD_BF_S(0, 1); ADD_BF_S(0, 2); ADD_BF_S(0, 3);
    ADD_BF_S(1, 0); ADD_BF_S(1, 1); ADD_BF_S(1, 2); ADD_BF_S(1, 3);
    ADD_BF_S(2, 0); ADD_BF_S(2, 1); ADD_BF_S(2, 2); ADD_BF_S(2, 3);
    ADD_BF_S(3, 0); ADD_BF_S(3, 1); ADD_BF_S(3, 2); ADD_BF_S(3, 3);

    ADD_LF(0);
    ADD_LF(1);
    ADD_LF(2);
    ADD_LF(3);

//    ADD_LF_S(1);
//    ADD_LF_S(2);

    // this is necessary, so that we can use Python from forms.cpp:
    if (import_hermes2d___hermes2d())
        throw std::runtime_error("hermes2d failed to import.");
}
