#include "numerical_flux.h"
#include "forms.h"

double matrix_R(int i, int j, double w0, double w1, double w3, double w4)
{
    double rho = w0;
    double u = w1/w0;
    double w = w3/w0;
    double E = w4;
    double v2 = u*u+w*w;
    double p = (kappa-1)*(E - rho*v2/2);
    double c = sqrt(kappa*p/rho);
    if (i == 0 && j == 0)
        return 1;
    else if (i == 0 && j == 1)
        return 1;
    else if (i == 0 && j == 2)
        return 1;
    else if (i == 0 && j == 3)
        return 1;

    else if (i == 1 && j == 0)
        return u-c;
    else if (i == 1 && j == 1)
        return u;
    else if (i == 1 && j == 2)
        return u;
    else if (i == 1 && j == 3)
        return u+c;

    else if (i == 2 && j == 0)
        return w;
    else if (i == 2 && j == 1)
        return w;
    else if (i == 2 && j == 2)
        return w-c;
    else if (i == 2 && j == 3)
        return w;

    else if (i == 3 && j == 0)
        return v2/2 + c*c/(kappa-1) - u*c;
    else if (i == 3 && j == 1)
        return v2/2;
    else if (i == 3 && j == 2)
        return v2/2 - w*c;
    else if (i == 3 && j == 3)
        return v2/2 + c*c/(kappa-1) + u*c;

    printf("i=%d, j=%d;\n", i, j);
    error("Invalid index.");
}

double matrix_R_inv(int i, int j, double w0, double w1, double w3, double w4)
{
    double rho = w0;
    double u = w1/w0;
    double w = w3/w0;
    double E = w4;
    double v2 = u*u+w*w;
    double p = (kappa-1)*(E - rho*v2/2);
    double c = sqrt(kappa*p/rho);
    double result;
    if (i == 0 && j == 0)
        result = ((kappa-1)*v2/2 + u*c)/2;
    else if (i == 0 && j == 1)
        result = -(c+u*(kappa-1))/2;
    else if (i == 0 && j == 2)
        result = -w*(kappa-1)/2;
    else if (i == 0 && j == 3)
        result = (kappa-1)/2;

    else if (i == 1 && j == 0)
        result = c*c-c*w-(kappa-1)*v2/2;
    else if (i == 1 && j == 1)
        result = u*(kappa-1);
    else if (i == 1 && j == 2)
        result = c+w*(kappa-1);
    else if (i == 1 && j == 3)
        result = 1-kappa;

    else if (i == 2 && j == 0)
        result = w*c;
    else if (i == 2 && j == 1)
        result = 0;
    else if (i == 2 && j == 2)
        result = -c;
    else if (i == 2 && j == 3)
        result = 0;

    else if (i == 3 && j == 0)
        result = ((kappa-1)*v2/2 - u*c)/2;
    else if (i == 3 && j == 1)
        result = (c-u*(kappa-1))/2;
    else if (i == 3 && j == 2)
        result = -w*(kappa-1)/2;
    else if (i == 3 && j == 3)
        result = (kappa-1)/2;
    else {
        printf("i=%d, j=%d;\n", i, j);
        error("Invalid index.");
    }
    return result/(c*c);
}

double matrix_D_minus(int i, int j, double w0, double w1, double w3, double w4)
{
    double rho = w0;
    double u = w1/w0;
    double w = w3/w0;
    double E = w4;
    double v2 = u*u+w*w;
    double p = (kappa-1)*(E - rho*v2/2);
    double c = sqrt(kappa*p/rho);
    double u_diag = 0;
    if (u < 0)
        u_diag = u;
    if (i == 0 && j == 0)
        return u-c;
    else if (i == 0 && j == 1)
        return 0;
    else if (i == 0 && j == 2)
        return 0;
    else if (i == 0 && j == 3)
        return 0;

    else if (i == 1 && j == 0)
        return 0;
    else if (i == 1 && j == 1)
        return u_diag;
    else if (i == 1 && j == 2)
        return 0;
    else if (i == 1 && j == 3)
        return 0;

    else if (i == 2 && j == 0)
        return 0;
    else if (i == 2 && j == 1)
        return 0;
    else if (i == 2 && j == 2)
        return u_diag;
    else if (i == 2 && j == 3)
        return 0;

    else if (i == 3 && j == 0)
        return 0;
    else if (i == 3 && j == 1)
        return 0;
    else if (i == 3 && j == 2)
        return 0;
    else if (i == 3 && j == 3)
        return 0;

    printf("i=%d, j=%d;\n", i, j);
    error("Invalid index.");
}

// multiplies two matrices
void dot(double result[4][4], double A[4][4], double B[4][4])
{
    for (int i=0; i < 4; i++)
        for (int j=0; j < 4; j++) {
            double sum=0;
            for (int k=0; k < 4; k++)
                sum += A[i][k] * B[k][j];
            result[i][j] = sum;
        }
}

// multiplies a matrix and a vector
void dot_vector(double result[4], double A[4][4], double B[4])
{
    for (int i=0; i < 4; i++) {
        double sum=0;
        for (int k=0; k < 4; k++)
            sum += A[i][k] * B[k];
        result[i] = sum;
    }
}

void T_rot(double result[4][4], double beta)
{
    for (int i=0; i < 4; i++)
        for (int j=0; j < 4; j++)
            result[i][j] = 0;
    result[0][0] = 1;
    result[1][1] = cos(beta);
    result[1][2] = sin(beta);
    result[2][1] = -sin(beta);
    result[2][2] = cos(beta);
    result[3][3] = 1;
}

void A_minus(double result[4][4], double w0, double w1, double w3, double w4)
{
    double _R[4][4];
    double _D_minus[4][4];
    double _R_inv[4][4];
    double _A_minus[4][4];
    double _tmp[4][4];
    for (int i=0; i < 4; i++)
        for (int j=0; j < 4; j++)
            _R[i][j] = matrix_R(i, j, w0, w1, w3, w4);
    for (int i=0; i < 4; i++)
        for (int j=0; j < 4; j++)
            _D_minus[i][j] = matrix_D_minus(i, j, w0, w1, w3, w4);
    for (int i=0; i < 4; i++)
        for (int j=0; j < 4; j++)
            _R_inv[i][j] = matrix_R_inv(i, j, w0, w1, w3, w4);
    dot(_tmp, _D_minus, _R_inv);
    dot(result, _R, _tmp);
}

void flux_riemann(double result[4], double w_l[4], double w_r[4])
{
    double _tmp1[4][4];
    double _tmp2[4][4];
    double _tmp3[4];
    double _tmp4[4];
    A_minus(_tmp1, w_r[0], w_r[1], w_r[2], w_r[3]);
    A_minus(_tmp2, w_l[0], w_l[1], w_l[2], w_l[3]);
    dot_vector(_tmp3, _tmp1, w_r);
    dot_vector(_tmp4, _tmp2, w_l);
    for (int i=0; i < 4; i++) {
        double _1 = f_x(i, w_l[0], w_l[1], w_l[2], w_l[3]);
        double _2 = _tmp3[i];
        double _3 = _tmp4[i];
        result[i] = _1 + _2 - _3;
    }
}

// calculates the inverted flux, for testing purposes
// it should return the same thing as flux_riemann(), only with minus sign
void flux_riemann_invert(double result[4], double w_l[4], double w_r[4])
{
    double m[4][4];
    double _w_l[4];
    double _w_r[4];
    double _tmp[4];
    T_rot(m, M_PI);
    dot_vector(_w_l, m, w_l);
    dot_vector(_w_r, m, w_r);
    flux_riemann(_tmp, _w_r, _w_l);
    T_rot(m, -M_PI);
    dot_vector(result, m, _tmp);
}
