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


Solution uprev;

#include "forms.h"

scalar residual(RealFunction* vi, RefMap* ri)
{
  RefMap* ru = uprev.get_refmap();
  return int_grad_u_grad_v(&uprev, vi, ru, ri) +
         int_u3(&uprev, ru);
}

scalar jacobian(RealFunction* vj, RealFunction* vi, RefMap* rj, RefMap* ri)
{
  RefMap* ru = uprev.get_refmap();
  return int_grad_u_grad_v(vj, vi, rj, ri) +
         3*int_u2_vj_vi(&uprev, vj, vi, ru, rj, ri); 
}


scalar bc_values(int marker, double x, double y)
{
  return 1.0 / hypot(x, y); 
}

scalar exact(double x, double y, scalar& dx, scalar& dy)
{
  return 1.0 / hypot(x, y);
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
  wf.add_biform(0, 0, jacobian, UNSYM, ANY, 1, &uprev);
  wf.add_liform(0, residual, ANY, 1, &uprev);
  
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);
  
  ScalarView view("Iteration", 0, 0, 850, 800);
  ScalarView errview("Error", 850, 0, 850, 800);
  
  uprev.set_zero(&mesh);
  
  do
  {
    Solution sln;
    nls.assemble();
    nls.solve(1, &sln);
    info("Residuum L2 norm: %g\n", nls.get_residuum_l2_norm());
    
    ExactSolution exsln(&mesh, exact);
    DiffFilter err(&sln, &exsln);
    view.show(&sln);
    errview.show(&err);
    errview.wait_for_keypress();
  }
  while (nls.get_residuum_l2_norm() > 1e-3);
  
  View::wait();
  return 0;
}
