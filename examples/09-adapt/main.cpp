#include "hermes2d.h"
#include "solver_umfpack.h"

// This example shows how to define different material parameters in
// various parts of the computational domain (see the definition of
// elements in the mesh file), and how to run automatic h- and
// hp-adaptivity with various parameters. The underlying problem is
// a planar model of an electrostatic micromotor. You may run the model
// with hp-adaptivity first (H_ONLY = false) and then with h-adaptivity
// (H_ONLY = true & P_INIT = 1, H_ONLY = true & P_INIT = 2). You can also
// check out the parameter ISO_ONLY that forbids/allows anisotropic
// element refinements (allowing anisotropic refinements is default
// in Hermes. Convergence graphs are saved (error estimate wrt. dof
// number and cpu time). Process the convergence files using Gnuplot,
// e.g., "gnuplot conv_dof.gp".
//
// PDE: -div[eps_r grad phi] = 0
//      eps_r = EPS1 in Omega_1 (surrounding air)
//      eps_r = EPS2 in Omega_2 (moving part of the motor)
//
// BC: phi = 0 V on Gamma_1 (left edge and also the rest of the outer
//               boundary
//     phi = VOLTAGE on Gamma_2 (boundary of stator)
//     
// The following parameters can be changed: 

const int P_INIT = 1;             // initial polynomial degree in mesh
const double THRESHOLD = 0.3;     // error threshold for element refinement
const int STRATEGY = 0;           // refinement strategy (0, 1, 2, 3 - see adapt_h1.cpp for explanation)
const bool H_ONLY = false;        // if H_ONLY == false then full hp-adaptivity takes place, otherwise
                                  // h-adaptivity is used. Use this parameter to check that indeed adaptive
                                  // hp-FEM converges much faster than adaptive h-FEM
const bool ISO_ONLY = false;      // when ISO_ONLY = true, only isotropic refinements are done,
                                  // otherwise also anisotropic refinements are allowed
const int MESH_REGULARITY = 1;    // specifies maximum allowed level of hanging nodes
                                  // -1 ... arbitrary level hangning nodes
                                  // 1, 2, 3,... means 1-irregular mesh, 2-irregular mesh, etc.
                                  // total regularization (0) is not supported in adaptivity
const double ERR_STOP = 0.1;      // stopping criterion for adaptivity (rel. error tolerance between the
                                  // reference and coarse solution in percent)
const int NDOF_STOP = 50000;      // adaptivity process stops when the number of degrees of freedom grows over
                                  // this limit. This is mainly to prevent h-adaptivity to go on forever.
const double EPS1 = 1.0;          // permittivity in Omega_1
const double EPS2 = 10.0;         // permittivity in Omega_2
const double VOLTAGE = 50.0;      // voltage on the stator


scalar bc_values(int marker, double x, double y)
{  
  return (marker == 2) ? VOLTAGE : 0.0;
}

scalar biform1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return EPS1 * int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar biform2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return EPS2 * int_grad_u_grad_v(fu, fv, ru, rv);
}


int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("motor.mesh");

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  H1Space space(&mesh, &shapeset);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
  
  WeakForm wf(1);
  wf.add_biform(0, 0, biform1, SYM, 1);
  wf.add_biform(0, 0, biform2, SYM, 2);
  
  ScalarView sview("Coarse solution", 0, 0, 600, 1000);
  VectorView gview("Gradient", 610, 0, 600, 1000);
  OrderView  oview("Polynomial orders", 1220, 0, 600, 1000);
  //gview.set_min_max_range(0.0, 400.0);
  
  Solution sln, rsln;
  UmfpackSolver solver;
  
  GnuplotGraph graph;
  graph.set_captions("Error Convergence for the Micromotor Problem", "Degrees of Freedom", "Error Estimate [%]");
  graph.add_row("error estimate", "-", "o");
  graph.set_log_y();

  GnuplotGraph graph_cpu;
  graph_cpu.set_captions("Error Convergence for the Micromotor Problem", "CPU Time", "Error Estimate [%]");
  graph_cpu.add_row("error estimate", "-", "o");
  graph_cpu.set_log_y();

  int it = 1;
  bool done = false;
  double cpu = 0.0;
  do
  {
    info("\n---- Iteration %d ---------------------------------------------\n", it++);
    begin_time();

    // enumerating basis functions
    space.assign_dofs();
  
    // solve the coarse problem
    LinSystem ls(&wf, &solver);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln);

    cpu += end_time();
    
    // view the solution -- this can be slow; for illustration only
    sview.show(&sln);
    gview.show(&sln, &sln, EPS_NORMAL, FN_DX_0, FN_DY_0);
    oview.show(&space);

    // solve the fine (reference) problem
    begin_time();
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &rsln);
    
    // calculate element errors and total error estimate
    H1OrthoHP hp(1, &space);
    double err_est = hp.calc_error(&sln, &rsln) * 100;
    info("Error estimate: %g%%", err_est);

    // adaptivity step
    if (err_est < ERR_STOP || ls.get_num_dofs() >= NDOF_STOP) done = true;
    else hp.adapt(THRESHOLD, STRATEGY, H_ONLY, ISO_ONLY, MESH_REGULARITY);
    cpu += end_time();

    // plotting convergence wrt. numer of dofs
    graph.add_values(0, space.get_num_dofs(), err_est);
    graph.save("conv_dof.gp");

    // plotting convergence wrt. cpu time
    graph_cpu.add_values(0, cpu, err_est);
    graph_cpu.save("conv_cpu.gp");
  }
  while (done == false);
  verbose("\nTotal running time: %g sec", cpu);
  
  // show the fine solution - this is the final result, about
  // one order of magnitude more accurate than the coarse solution
  sview.set_title("Fine solution");
  sview.show(&rsln); 
  gview.show(&rsln, &rsln, EPS_HIGH, FN_DX_0, FN_DY_0);

  View::wait();
  return 0;
}
