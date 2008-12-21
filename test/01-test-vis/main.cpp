// *** visualization test ***

#include "hermes2d.h"


const int order = 2;


int xvel_bc_type(int marker)
  { return (marker != 11) ? BC_ESSENTIAL : BC_NONE; }

scalar xvel_bc_value(int marker, double x, double y)
  { return (marker > 3) ? 1 : 0; }

int yvel_bc_type(int marker)
  { return (marker != 11) ? BC_ESSENTIAL : BC_NONE; }

int press_bc_type(int marker)
  { return BC_NONE; }

  

int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);
  
  Mesh mesh;
  mesh.load("airfoil.mesh");
  
  MeshView mv;
  mv.show(&mesh);
  
  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);
    
  // H1 spaces for velocities and pressure
  H1Space xvel(&mesh, &shapeset);
  H1Space yvel(&mesh, &shapeset);
  H1Space press(&mesh, &shapeset);

  // set velocity and pressure polynomial degrees
  xvel.set_uniform_order(order);
  yvel.set_uniform_order(order);
  press.set_uniform_order(order-1);

  // initialize boundary conditions
  xvel.set_bc_types(xvel_bc_type);
  xvel.set_bc_values(xvel_bc_value);
  yvel.set_bc_types(yvel_bc_type);
  press.set_bc_types(press_bc_type);

  // assign degrees of freedom
  int ndofs = 0;
  ndofs += xvel.assign_dofs(ndofs);
  ndofs += yvel.assign_dofs(ndofs);
  ndofs += press.assign_dofs(ndofs);
  printf("ndofs = %d\n", ndofs);

  // load solution vector from a file
  scalar* vec = new scalar[ndofs+1];
  FILE* f = fopen("sln200.dat", "rb");
  fread(vec, sizeof(double), ndofs+1, f);
  fclose(f);
  
  Solution xsln, ysln;
  xsln.set_space_and_pss(&xvel, &pss);
  ysln.set_space_and_pss(&yvel, &pss);
  xsln.set_solution_vector(vec, true);  
  ysln.set_solution_vector(vec, false);

  ScalarView view1("Scalar1", 200, 150, 1100, 800);
  view1.show(&xsln);
  ScalarView view2("Scalar2", 200, 150, 1100, 800);
  view2.show(&ysln);
  VectorView vview("Vector", 200, 150, 1100, 800);
  vview.show(&xsln, &ysln);
    
  hermes2d_finalize();
  return 0;
}
