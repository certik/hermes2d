//
//  Multi-mesh assembling test:
//
//  Domain:  (-1;1)^2,  Neumann: x=1, Dirichlet: elsewhere
//
//  Exact solution:  u0 = x^2 + y^2
//                   u1 = x^3
//
//  PDEs:  -\Delta u0 + u1 = f0
//          u0 - \Delta u1 = f1
//
//          f0 = 4 + x^3
//          f1 = x^2 + y^2 + 6x
//

#include "hermes2d.h"
#include "solver_umfpack.h"


int bc_types(int marker)
  { return (marker == 2) ? BC_NATURAL : BC_ESSENTIAL; }

double bc_values_0(int marker, double x, double y)
  { return x*x + y*y; }

double bc_values_1(int marker, double x, double y)
  { return x*x*x; }


scalar bilinear_form_unsym_0_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv); }

scalar bilinear_form_unsym_0_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_u_v(fu, fv, ru, rv); }

scalar bilinear_form_unsym_1_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_u_v(fu, fv, ru, rv); }

scalar bilinear_form_unsym_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv); }

  
double f_0(double x, double y)
  { return -4 + x*x*x; }
  
double f_1(double x, double y)
  { return x*x + y*y - 6*x; }

scalar linear_form_0(RealFunction* fv, RefMap* rv)
  { return int_F_v(f_0, fv, rv); }

scalar linear_form_1(RealFunction* fv, RefMap* rv)
  { return int_F_v(f_1, fv, rv); }

scalar linear_form_surf_0(RealFunction* fv, RefMap* rv, EdgePos* ep)
  { return (ep->marker == 2) ? 2.0 * surf_int_v(fv, rv, ep) : 0.0; }
  
scalar linear_form_surf_1(RealFunction* fv, RefMap* rv, EdgePos* ep)
  { return (ep->marker == 2) ? 3.0 * surf_int_v(fv, rv, ep) : 0.0; }

  
  
int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);
  
  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss0(&shapeset);
  PrecalcShapeset pss1(&shapeset);
  
  Mesh mesh0, mesh1;
  mesh0.load("square1.mesh");
  mesh1.copy(&mesh0);
  
  /*mesh0.refine_all_elements();
  mesh1.refine_all_elements();
  mesh0.refine_all_elements();
  mesh1.refine_all_elements();*/
  
  mesh1.refine_towards_vertex(3, 5);

  H1Space space0(&mesh0, &shapeset);
  space0.set_bc_types(bc_types);
  space0.set_bc_values(bc_values_0);
  space0.set_uniform_order(2);
  int ndofs = space0.assign_dofs();

  H1Space space1(&mesh1, &shapeset);
  space1.set_bc_types(bc_types);
  space1.set_bc_values(bc_values_1);
  space1.set_uniform_order(3);
  ndofs += space1.assign_dofs(ndofs);

  WeakForm wf(2);
  wf.add_biform(0, 0, bilinear_form_unsym_0_0);
  wf.add_biform(0, 1, bilinear_form_unsym_0_1);
  wf.add_biform(1, 0, bilinear_form_unsym_1_0);
  wf.add_biform(1, 1, bilinear_form_unsym_1_1);
  wf.add_liform(0, linear_form_0);
  wf.add_liform(1, linear_form_1);
  wf.add_liform_surf(0, linear_form_surf_0);
  wf.add_liform_surf(1, linear_form_surf_1);
  
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(2, &space0, &space1);
  sys.set_pss(2, &pss0, &pss1);
  
  Solution sln0, sln1;
  sys.assemble();
  sys.solve(2, &sln0, &sln1);

  ScalarView view0("Component 0", 10, 200, 780, 780);
  ScalarView view1("Component 1", 800, 200, 780, 780);
  
  view0.show(&sln0);
  view1.show(&sln1);

  hermes2d_finalize();
  return 0;
}

