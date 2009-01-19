// 
//  Nonlinear solver test:
//
//  PDE:  -\Delta u + u^3 = 0
//
//  Domain: annulus centered at origin, with inner radius A = 0.05 and
//  outer radius B = 0.5.
//
//  BC:  u = 1/A  on inner boundary
//       du/dn = -1/B^2  on outer boundary
//
//  Exact solution:  u = 1 / sqrt(x^2 + y^2)
//

#include "hermes2d.h"
#include "solver_umfpack.h"

scalar bc_values(int marker, double x, double y)
{
  return 1 / hypot(x, y); 
}


#include "forms.h"

Solution uprev;

scalar residual(RealFunction* vi, RefMap* ri)
{
  RefMap* ru = uprev.get_refmap();
  return int_grad_u_grad_v(uprev, vi, ru, ri) +
         int_u3(uprev, ru);
}

scalar jacobian(RealFunction* vj, RealFunction* vi, RefMap* rj, RefMap* ri)
{
  RefMap* ru = uprev.get_refmap();
  return int_grad_u_grad_v(fu, fv, ru, rv) +
         3*int_u2_vj_vi(uprev, vj, vi, ru, rj, ri); 
}


scalar biform(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}


scalar exact(double x, double y, scalar& dx, scalar& dy)
{
  return 1 / hypot(x, y);
}


int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("annulus.mesh");
  mesh.refine_all_elements(1);
  mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, 3);

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  H1Space space(&mesh, &shapeset);
  space.set_bc_values(bc_values);
  space.set_uniform_order(3);
  space.assign_dofs();
  
  WeakForm wf(1);
  wf.add_biform(0, 0, biform);
  
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);
  
  Solution sln;
  sys.assemble();
  sys.solve(1, &sln);
  
  ScalarView view;
  view.show(&sln);
  
  ScalarView view2;
  ExactSolution esln(&mesh, exact);
  view2.show(&esln);
  
  View::wait();
  return 0;
}

