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
    w0_init_num = rho_z(0)/rho_r;
    return w0_init_num;
}

scalar w1_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
    w1_init_num = rho_z(0)/rho_r * (20/u_r);
    return w1_init_num;
}

scalar w3_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
    w3_init_num = rho_z(0)/rho_r * (0/u_r);
    return w3_init_num;
}


scalar w4_init(double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
    w4_init_num = rho_z(0) * T_z(0) * c_v / E_r;
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
/*
    double w_l[4];
    // XXX: set w_l and w_r to the prev solution
    double w_r[4];
    // set the optional BC into w_r
    if (e->marker == marker_top || e->marker == marker_bottom) {
        // the z-velocity is 0:
        w_r[2] = 0;
    }
    if (e->marker == marker_left || e->marker == marker_right) {
        // the x-velocity is 1 m/s:
        w_r[1] = 1./u_r;
    }

    double alpha = 0; // XXX rotate things appropriately

    double m[4][4];
    double _w_l[4];
    double _w_r[4];
    double _tmp[4];
    double flux[4];
    T_rot(m, alpha);
    dot_vector(_w_l, m, w_l);
    dot_vector(_w_r, m, w_r);
    flux_riemann(_tmp, _w_l, _w_r);
    T_rot(m, -alpha);
    dot_vector(flux, m, _tmp);
*/

    double w0, w1, w3, w4;

    Scalar result = 0;
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
        //printf("c = %f; M = %f\n", _c, );
        if (e->marker == marker_top || e->marker == marker_bottom) {
            double un = _u*e->nx[i] + _w*e->ny[i];
            //printf("normal part: %f\n", un);
            //printf("BC: %f %f \n", w1, w3);
            _u = _u - 2 * un * e->nx[i];
            _w = _w - 2 * un * e->ny[i];
            w1 = _u * w0;
            w3 = _w * w0;
            /*
            w1 = -0.5*10;
            if (e->marker == marker_top)
                w3 = -0.5*10;
            else
                w3 = 0.5*10;
            */
            //printf("BC(%d): %f %f %f %f \n", e->marker, w1, w3, e->nx[i], e->ny[i]);
        }
        else {
            double un = _u*e->nx[i] + _w*e->ny[i];
            if (un > 0) {
                // outlet
                if (_M >= 1.) {
                    // supersonic
                    // we take everything from inside
                }
                else {
                    // subsonic
                    // take p from the outside state
                    w4 = p_init_num * c_v / R + (w1*w1+w3*w3)/(2*w0);
                    //w4 = rho_z(0) * T_z(0) * c_v / E_r;
                }
            }
            else {
                // inlet
                if (_M >= 1.) {
                    // supersonic
                    // take everything from outside
                    w0 = w0_init_num;
                    w1 = w1_init_num;
                    w3 = w3_init_num;
                    w4 = w4_init_num;
                }
                else {
                    // subsonic
                    // take rho, u, w from the outside state
                    w0 = w0_init_num;
                    w1 = w1_init_num;
                    w3 = w3_init_num;
                }
                /*
                printf("left: %f %f %f, %f %f %f\n", w0, w1, w3, rho, u*rho,
                        w*rho);
                        */
            }
        }
        result += wt[i] * (
                A_x(_i, _j, w0, w1, w3, w4) * e->nx[i]
                +
                A_z(_i, _j, w0, w1, w3, w4) * e->ny[i]
                ) * u->val[i] * v->val[i];
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

#define ADD_BF(i, j) wf.add_biform(i, j, callback_bf(B_##i##j), UNSYM, ANY, 4, &w0_prev, &w1_prev, &w3_prev, &w4_prev)

#define ADD_BF_S(i, j) wf.add_biform_surf(i, j, callback_bf(S_##i##j), ANY, 4, &w0_prev, &w1_prev, &w3_prev, &w4_prev)

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

    wf.add_liform(0, callback(l_0), ANY, 4, &w0_prev, &w1_prev, &w3_prev,
            &w4_prev);
    wf.add_liform(1, callback(l_1), ANY, 4, &w0_prev, &w1_prev, &w3_prev,
            &w4_prev);
    wf.add_liform(2, callback(l_2), ANY, 4, &w0_prev, &w1_prev, &w3_prev,
            &w4_prev);
    wf.add_liform(3, callback(l_3), ANY, 4, &w0_prev, &w1_prev, &w3_prev,
            &w4_prev);
}
