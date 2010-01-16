#ifndef __HERMES2D_SEA_BREEZE_NUMERICAL_FLUX_H
#define __HERMES2D_SEA_BREEZE_NUMERICAL_FLUX_H

#include "hermes2d.h"
#include "params.h"

double matrix_R(int i, int j, double w0, double w1, double w3, double w4);
double matrix_R_inv(int i, int j, double w0, double w1, double w3, double w4);
double matrix_D_minus(int i, int j, double w0, double w1, double w3, double w4);
void flux_riemann(double result[4], double w_l[4], double w_r[4]);

#endif
