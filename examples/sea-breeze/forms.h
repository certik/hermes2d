#ifndef __HERMES2D_SEA_BREEZE_FORMS_H
#define __HERMES2D_SEA_BREEZE_FORMS_H

#include "hermes2d.h"

void register_bc(H1Space &s0, H1Space &s1, H1Space &s3, H1Space &s4);
void set_ic(Mesh &mesh, Solution &w0, Solution &w1, Solution &w3, Solution &w4);
void register_forms(WeakForm &wf, Solution &w0_prev, Solution &w1_prev,
        Solution &w3_prev, Solution &w4_prev);

extern const double TAU;

#endif
