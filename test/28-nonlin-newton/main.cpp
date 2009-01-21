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

const double A = 0.05;
const double B = 0.5;


Solution uprev;

#include "forms.h"

scalar residual(RealFunction* vi, RefMap* ri)
{
  return int_grad_U_grad_v(&uprev, vi, ri, ri) +
         int_u3_v(&uprev, vi, ri, ri);
}

scalar residual_surf(RealFunction* vi, RefMap* ri, EdgePos* ep)
{
  return 1.0/(sqr(B)) * surf_int_v(vi, ri, ep);
}


scalar jacobian(RealFunction* vj, RealFunction* vi, RefMap* rj, RefMap* ri)
{
  RefMap* ru = uprev.get_refmap();
  return int_grad_u_grad_v(vj, vi, rj, ri) +
         3*int_u2_vj_vi(&uprev, vj, vi, ru, rj, ri);
}


int bc_types(int marker)
{
  if (marker == 1)
    return BC_ESSENTIAL;
  return BC_NATURAL;
}

scalar bc_values(int marker, double x, double y)
{
  if (marker == 1) return 1.0 / A; 
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
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(3);
  space.assign_dofs();
  
  WeakForm wf(1);
  wf.add_biform(0, 0, jacobian, UNSYM, ANY, 1, &uprev);
  wf.add_liform(0, residual, ANY, 1, &uprev);
  wf.add_liform_surf(0, residual_surf, 2);
  
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);
  
  ScalarView view("Iteration", 0, 0, 880, 800);
  ScalarView errview("Error", 900, 0, 880, 800);
  
  // uprev must contain the dirichlet lift at the beginning
  // TODO: make a new function of Solution for this
  scalar vec[space.get_num_dofs()];
  memset(vec, 0, sizeof(vec));
  uprev.set_fe_solution(&space, &pss, vec);
  
  int it = 1;
  double l2e;
  do
  {
    info("\n---- Iteration %d ---------------------------------------\n", it++);
    
    Solution sln;
    nls.assemble();
    nls.solve(1, &sln);
    info("Residuum L2 norm: %g\n", nls.get_residuum_l2_norm());
    
    view.show(&sln);
    
    ExactSolution exsln(&mesh, exact);
    DiffFilter err(&sln, &exsln);
    errview.show(&err);
    errview.wait_for_keypress();
    
    l2e = l2_error(&sln, &exsln);
    info("L2 error: %g", l2e);
    
    uprev = sln;
  }
  while (nls.get_residuum_l2_norm() > 1e-4);
  //while (l2e > 1e-5);
  
  View::wait();
  return 0;
}
