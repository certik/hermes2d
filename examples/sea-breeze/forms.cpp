#include "forms.h"
#include "numerical_flux.h"

#include "python_api.h"

const double TAU = 1e-5;  // this is in seconds

const double T_0 = T_z(0);
const double p_0 = p_z(0);
const double rho_0 = rho_z(0);

BCType s0_bc_type(int marker) {
    return BC_NATURAL;
}

BCType s1_bc_type(int marker) {
    return BC_NATURAL;
}

BCType s3_bc_type(int marker) {
    return BC_NATURAL;
}

BCType s4_bc_type(int marker) {
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
    w1_init_num = 50;
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
    w4_init_num = 1e5;
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

double f_x_via_A(int i, double w0, double w1, double w3, double w4)
{
    double r=0;
    for (int j=0; j<4; j++) {
        r += A_x(i, j, w0, w1, w3, w4) * w(j, w0, w1, w3, w4);
    }
    return r;
}

double f_z_via_A(int i, double w0, double w1, double w3, double w4)
{
    double r=0;
    for (int j=0; j<4; j++) {
        r += A_z(i, j, w0, w1, w3, w4) * w(j, w0, w1, w3, w4);
    }
    return r;
}

double test1(double w0, double w1, double w3, double w4)
{
    for (int i=0; i<4; i++) {
        double r=f_x_via_A(i, w0, w1, w3, w4);
        r -= f_x(i, w0, w1, w3, w4);
        printf("result: %d: %f\n", i, r);
    }
}

double test3(double w0, double w1, double w3, double w4)
{
    for (int i=0; i<4; i++) {
        double r=f_z_via_A(i, w0, w1, w3, w4);
        r -= f_z(i, w0, w1, w3, w4);
        printf("result: %d: %f\n", i, r);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////

int counter=0;

template<typename Real, typename Scalar>
Scalar B_ij(int _i, int _j, int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    assert(u->val != NULL);
    assert(v->dx != NULL);
    assert(v->dy != NULL);
    counter++;
    /*
    insert_object("fn", array_double_c2numpy(ext->fn[3]->val, n));
    insert_object("v", array_double_c2numpy(v->val, n));
    insert_object("x", array_double_c2numpy(e->x, n));
    insert_object("y", array_double_c2numpy(e->y, n));
    insert_object("tau", double_c2py(TAU));
    */
    Python p;
    p.push("i", c2py_int(_i));
    p.push("j", c2py_int(_j));
    p.push("counter", c2py_int(counter));
    //p.exec("print '-'*40");
    //p.exec("print 'B_ij called, counter=', counter");
    //p.exec("print 'i=%d j=%d' % (i, j)");
    double delta_ij;
    if (_i == _j)
        delta_ij = 1;
    else
        delta_ij = 0;
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (
                (u->val[i] * v->val[i]) * delta_ij / TAU /* +
                -A_x(_i, _j, ext->fn[0]->val[i],
                    ext->fn[1]->val[i],
                    ext->fn[2]->val[i],
                    ext->fn[3]->val[i])
                * u->val[i] * v->dx[i] +
                -A_z(_i, _j, ext->fn[0]->val[i],
                    ext->fn[1]->val[i],
                    ext->fn[2]->val[i],
                    ext->fn[3]->val[i])
                * u->val[i] * v->dy[i] */
                );
    return result;
}

template<typename Real, typename Scalar>
Scalar S_ij(int _i, int _j, int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    double w0, w1, w3, w4;
    w0 = ext->fn[0]->val[0];
    w1 = ext->fn[1]->val[0];
    w3 = ext->fn[2]->val[0];
    w4 = ext->fn[3]->val[0];

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
    return result;
}

template<typename Real, typename Scalar>
Scalar s_i(int _i, int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    //printf("n=%d\n", n);
    assert(n < 20);
    double w0, w1, w3, w4;
    Scalar result = 0;
    //printf("BC: n=%d; marker=%d\n", n, e->marker);
    for (int i = 0; i < n; i++) {
        double w_l[4];
        double w_r[4];
        double flux[4];
        double nx=e->nx[i], ny=e->ny[i];

        // this has to be fixed in hermes to return the correct w_l, w_r:
        for (int j=0; j < 4; j++)
            w_l[j] = ext->fn[j]->val[i];
        if (e->marker == marker_top || e->marker == marker_bottom) {
            double w0, w1, w3, w4;
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

                double alpha = atan2(ny, nx);
                double mat_rot_inv[4][4];
                double flux_local[4];
                flux_local[0] = 0;
                flux_local[1] = _p;
                flux_local[2] = 0;
                flux_local[3] = 0;
                T_rot(mat_rot_inv, -alpha);
                dot_vector(flux, mat_rot_inv, flux_local);
        } else if (e->marker == marker_left || e->marker == marker_right) {
            w_r[0] = w0_init_num;
            w_r[1] = w1_init_num;
            w_r[2] = w3_init_num;
            w_r[3] = w4_init_num;
            numerical_flux(flux, w_l, w_r, nx, ny);
        } else {
            // get from the neighbor element
            for (int j=0; j < 4; j++)
                w_r[j] = ext->fn2[j]->val[i];
            //printf("w_r: %f %f %f %f\n", w_r[0], w_r[1], w_r[2], w_r[3]);
            /*
            w_r[0] = w0_init_num;
            w_r[1] = w1_init_num;
            w_r[2] = w3_init_num;
            w_r[3] = w4_init_num;
            */
            numerical_flux(flux, w_l, w_r, nx, ny);
        }

        result += wt[i] * flux[_i] * v->val[i];
    }
    return result;
}

template<typename Real, typename Scalar>
Scalar l_i(int _i, int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    double w0, w1, w3, w4;
    Scalar result = 0;
    Scalar result_flux = 0;
    for (int i = 0; i < n; i++) {
        w0 = ext->fn[0]->val[i];
        w1 = ext->fn[1]->val[i];
        w3 = ext->fn[2]->val[i];
        w4 = ext->fn[3]->val[i];
        result += wt[i] * (
                    + ext->fn[_i]->val[i] / TAU * v->val[i]
                );
        result_flux += wt[i] * (
                      f_x_via_A(_i, w0, w1, w3, w4) * v->dx[i]
                    + f_z_via_A(_i, w0, w1, w3, w4) * v->dy[i]
                );
    }
    //printf("result=%f, result_flux=%f\n", result, result_flux);
    return result+result_flux;
    /*
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (ext->fn[0]->val[i] * v->val[i]) / TAU;
    return result;
    */
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
    return l_i(0, n, wt, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar l_1(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return l_i(1, n, wt, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar l_2(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return l_i(2, n, wt, v, e, ext);
}

//int counter=0;

template<typename Real, typename Scalar>
Scalar l_3(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return l_i(3, n, wt, v, e, ext);
    /*
    counter++;
    insert_object("fn", array_double_c2numpy(ext->fn[3]->val, n));
    insert_object("v", array_double_c2numpy(v->val, n));
    insert_object("x", array_double_c2numpy(e->x, n));
    insert_object("y", array_double_c2numpy(e->y, n));
    insert_object("tau", double_c2py(TAU));
    insert_object("counter", int_c2py(counter));
    cmd("print '-'*40");
    cmd("print 'l_3 called, counter=', counter");
    */
    //cmd("print 'fn', fn");
    //cmd("print 'v', v");
    //cmd("print 'x', x");
    //cmd("print 'y', y");
    //cmd("import util");
    //cmd("util.plotxy(x, y, v, counter)");
    //cmd("print tau");
    /*
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (ext->fn[3]->val[i] * v->val[i]) / TAU;
        */
    /*
    insert_object("result", double_c2py(result));
    cmd("print result");
    cmd("print '-'*40");
    */
    //return result;
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
Scalar s_0(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return s_i(0, n, wt, v, e, ext);
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
Scalar s_3(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return s_i(3, n, wt, v, e, ext);
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

#define callback_lf(a) a<double, scalar>, L_order

#define ADD_BF(i, j) wf.add_biform(i, j, callback_bf(B_##i##j), H2D_UNSYM, H2D_ANY, 4, &w0_prev, &w1_prev, &w3_prev, &w4_prev)

#define ADD_BF_S(i, j) wf.add_biform_surf(i, j, callback_bf(S_##i##j), H2D_ANY_EDGE, 4, &w0_prev, &w1_prev, &w3_prev, &w4_prev)

#define ADD_LF(i) wf.add_liform(i, callback_lf(l_##i), H2D_ANY, 4, &w0_prev, &w1_prev, &w3_prev, &w4_prev);

#define ADD_LF_S(i) wf.add_liform_surf(i, callback_lf_s(s_##i), H2D_ANY_EDGE, 4, &w0_prev, &w1_prev, &w3_prev, &w4_prev);

void register_forms(WeakForm &wf, Solution &w0_prev, Solution &w1_prev,
        Solution &w3_prev, Solution &w4_prev)
{
    ADD_BF(0, 0); ADD_BF(0, 1); ADD_BF(0, 2); ADD_BF(0, 3);
    ADD_BF(1, 0); ADD_BF(1, 1); ADD_BF(1, 2); ADD_BF(1, 3);
    ADD_BF(2, 0); ADD_BF(2, 1); ADD_BF(2, 2); ADD_BF(2, 3);
    ADD_BF(3, 0); ADD_BF(3, 1); ADD_BF(3, 2); ADD_BF(3, 3);

    /*
    ADD_BF_S(0, 0); ADD_BF_S(0, 1); ADD_BF_S(0, 2); ADD_BF_S(0, 3);
    ADD_BF_S(1, 0); ADD_BF_S(1, 1); ADD_BF_S(1, 2); ADD_BF_S(1, 3);
    ADD_BF_S(2, 0); ADD_BF_S(2, 1); ADD_BF_S(2, 2); ADD_BF_S(2, 3);
    ADD_BF_S(3, 0); ADD_BF_S(3, 1); ADD_BF_S(3, 2); ADD_BF_S(3, 3);
    */

    ADD_LF(0);
    ADD_LF(1);
    ADD_LF(2);
    ADD_LF(3);

    ADD_LF_S(0);
    ADD_LF_S(1);
    ADD_LF_S(2);
    ADD_LF_S(3);
}
