//  Newton BC example:
//
//  (Stationary heat transfer with constant temperature on Gamma_3, insulated
//   walls Gamma_2, Gamma_4 and Netwon-type cooling on Gamma_1).
//
//  PDE:  -\Delta u = 0
//
//  BC:   u = t1              on  Gamma_3
//        du/dn = 0           on  Gamma_2 & Gamma_4
//        du/dn = h*(u - t0)  on  Gamma_1
//
//  (Note that the last BC can be written in the form  du/dn - h*u = -h*t0 )

#include "hermes2d.h"
#include "solver_umfpack.h"

const double t1 = 100.0;
const double t0 = 20.0;
const double h  = 1.0;


int bc_types(int marker)
  { return (marker == 3) ? BC_ESSENTIAL : BC_NATURAL; }

scalar bc_values(int marker, double x, double y)
  { return (marker == 3) ? t1 : 0.0; }


scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv); }

scalar bilinear_form_surf_Gamma_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv, EdgePos* ep)
  { return h * surf_int_u_v(fu, fv, ru, rv, ep); }

scalar linear_form_surf_Gamma_1(RealFunction* fv, RefMap* rv, EdgePos* ep)
  { return t0 * h * surf_int_v(fv, rv, ep); }


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(6);
  space.assign_dofs();

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form);
  wf.add_biform_surf(0, 0, bilinear_form_surf_Gamma_1, 1);
  wf.add_liform_surf(0, linear_form_surf_Gamma_1, 1);
  
  // initialize the linear system and solver
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);
  
  // assemble the stiffness matrix and solve the system
  Solution sln;
  sys.assemble();
  sys.solve(1, &sln);
  
  // visualize the solution
  ScalarView view("Solution");
  view.show(&sln);
  
  View::wait();
  return 0;
}
