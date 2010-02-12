#include <cmath>
#include <hermes2d.h>
#include "functions.h"

//x2
scalar func_x2_val(double x, double y) { return x*x; };
scalar func_x2_dx(double x, double y) { return 2*x; };
scalar func_x2_dy(double x, double y) { return 0; };

//abs(y)
scalar func_absy_val(double x, double y) { return std::abs(y); };
scalar func_absy_dx(double x, double y) { return 0; };
scalar func_absy_dy(double x, double y) { if (y < 0) return -1; else if (y > 0) return 1; else return 0; };

//abs(x)*abs(y)
scalar func_absx_absy_val(double x, double y) { return std::abs(x)*std::abs(x); };
scalar func_absx_absy_dx(double x, double y) { if (x < 0) return -std::abs(y); else if (x > 0) return std::abs(y); else return 0; };
scalar func_absx_absy_dy(double x, double y) { if (y < 0) return -std::abs(x); else if (y > 0) return std::abs(x); else return 0; };

//x2y2
scalar func_x2y2_val(double x, double y) { return x*x * y*y; };
scalar func_x2y2_dx(double x, double y) { return 2*x * y*y; };
scalar func_x2y2_dy(double x, double y) { return x*x * 2*y; };

//x3y1
scalar func_x3y1_val(double x, double y) { return x*x*x * y; };
scalar func_x3y1_dx(double x, double y) { return 3*x*x * y; };
scalar func_x3y1_dy(double x, double y) { return x*x*x; };

//x3y4
scalar func_x3y4_val(double x, double y) { return x*x*x * y*y*y*y; };
scalar func_x3y4_dx(double x, double y) { return 3*x*x * y*y*y*y; };
scalar func_x3y4_dy(double x, double y) { return x*x*x * 4*y*y*y; };

