#include "hermes2d.h"
#include "solver_umfpack.h"

const double epsilon1 = 1.0;  // permittivity in Omega_1
const double epsilon2 = 10.0; // permittivity in Omega_2

const double tol = 0.01; // error tolerance in percent
const double thr = 0.3;  // error threshold for element refinement


scalar bc_values(int marker, double x, double y)
{  
  return (marker == 2) ? 50.0 : 0.0;
}

scalar biform1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return epsilon1 * int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar biform2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return epsilon2 * int_grad_u_grad_v(fu, fv, ru, rv);
}


int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("motor.mesh");

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  H1Space space(&mesh, &shapeset);
  space.set_bc_values(bc_values);
  space.set_uniform_order(1);
  
  WeakForm wf(1);
  wf.add_biform(0, 0, biform1, SYM, 1);
  wf.add_biform(0, 0, biform2, SYM, 2);
  
  ScalarView sview("Coarse solution", 0, 0, 600, 1000);
  VectorView gview("Gradient", 610, 0, 600, 1000);
  OrderView  oview("Polynomial orders", 1220, 0, 600, 1000);
  gview.set_min_max_range(0.0, 400.0);
  
  Solution sln, rsln;
  UmfpackSolver solver;
  
  GnuplotGraph graph("Error convergence");
  graph.add_row("hp-adaptivity", "k", "-","o");
  graph.set_log_y();

  int it = 1;
  while (1)
  {
    info("\n---- Iteration %d ---------------------------------------------\n", it++);
    space.assign_dofs();
  
    // solve the coarse problem
    LinSystem ls(&wf, &solver);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln);
    
    // view the solution -- this can be slow; for illustration only
    sview.show(&sln);
    gview.show(&sln, &sln, EPS_NORMAL, FN_DX_0, FN_DY_0);
    oview.show(&space);
    info("");

    // solve the fine (reference) problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &rsln);
    
    // calculate errors and adapt the solution
    H1OrthoHP hp(1, &space);
    double error = hp.calc_error(&sln, &rsln) * 100;
    graph.add_values(0, space.get_num_dofs(), error);
    graph.save("motor_conv_hp.txt");
    if (error < tol) break;
    hp.adapt(thr);
  }
  
  // show the fine solution - this is the final result, about
  // one order of magnitude more accurate than the coarse solution
  sview.set_title("Fine solution");
  sview.show(&rsln); 
  gview.show(&rsln, &rsln, EPS_HIGH, FN_DX_0, FN_DY_0);
  
  View::wait();
  return 0;
}
