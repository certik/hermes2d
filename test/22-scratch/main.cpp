#include "hermes2d.h"
#include "solver_umfpack.h"


scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return int_v(fv, rv);
}


int main(int argc, char* argv[])
{
  /*Mesh mesh;
  mesh.load_new("domain.mesh");*/
  
  /*mesh.load("square.mesh");
  mesh.refine_all_elements();
  mesh.refine_element(1);
  mesh.refine_all_elements();
  mesh.save("saved.mesh");*/
  
  //mesh.load("/home/jakub/hermes2d/examples/01-mesh/domain.mesh");
  //mesh.save("domain.mesh");
    
  /*MeshView mv;
  mv.show(&mesh);*/
  
  
  /*mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();*/
  
  /*H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  H1Space space(&mesh, &shapeset);
  space.set_uniform_order(10);
  space.assign_dofs();
  
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form, 1);
  wf.add_liform(0, linear_form);
  
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);
  
  Solution sln;
  sys.assemble();
  sys.solve(1, &sln);
  
  ScalarView view1;
  //view1.show_contours(0.02);
  //view1.show(&sln);
  view1.show(&sln, EPS_HIGH, FN_VAL_0);*/
  
  Solution sln;
  sln.load("bessel.sln.gz");
  
  Solution cpy;
  cpy.copy(&sln);
  
  RealFilter real(&cpy);
  MagFilter mag(&real);
  ScalarView view3("Magnitude of real(E)", 100, 50, 1000, 900);
  view3.set_min_max_range(0, 1);
  view3.show(&mag);
  
  View::wait();
  return 0;
}

