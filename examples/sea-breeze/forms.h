#ifndef __HERMES2D_SEA_BREEZE_FORMS_H
#define __HERMES2D_SEA_BREEZE_FORMS_H

#include "hermes2d.h"
#include "params.h"

//  boundary markers
#define marker_bottom 1
#define marker_right 2
#define marker_top 3
#define marker_left 4

void set_iteration(int i);
void register_bc(H1Space &s0, H1Space &s1, H1Space &s3, H1Space &s4);
void set_ic(Mesh &mesh, Solution &w0, Solution &w1, Solution &w3, Solution &w4);
void register_forms(WeakForm &wf, Solution &w0_prev, Solution &w1_prev,
        Solution &w3_prev, Solution &w4_prev);

double f_x(int i, double w0, double w1, double w3, double w4);
double f_z(int i, double w0, double w1, double w3, double w4);

double A_x(int i, int j, double w0, double w1, double w3, double w4);
double A_z(int i, int j, double w0, double w1, double w3, double w4);

extern const double TAU;

#endif
