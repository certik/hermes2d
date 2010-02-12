#ifndef __H2D_TEST_CAND_PROJ_H
#define __H2D_TEST_CAND_PROJ_H

typedef scalar (*ValueFunction) (double x, double y); ///< A function that returns a value.

//x2
extern scalar func_x2_val(double x, double y);
extern scalar func_x2_dx(double x, double y);
extern scalar func_x2_dy(double x, double y);

//abs(y)
extern scalar func_absy_val(double x, double y);
extern scalar func_absy_dx(double x, double y);
extern scalar func_absy_dy(double x, double y);

//abs(x)*abs(y)
extern scalar func_absx_absy_val(double x, double y);
extern scalar func_absx_absy_dx(double x, double y);
extern scalar func_absx_absy_dy(double x, double y);

//x2y2
extern scalar func_x2y2_val(double x, double y);
extern scalar func_x2y2_dx(double x, double y);
extern scalar func_x2y2_dy(double x, double y);

//x3y1
extern scalar func_x3y1_val(double x, double y);
extern scalar func_x3y1_dx(double x, double y);
extern scalar func_x3y1_dy(double x, double y);

//x3y4
extern scalar func_x3y4_val(double x, double y);
extern scalar func_x3y4_dx(double x, double y);
extern scalar func_x3y4_dy(double x, double y);

#endif

