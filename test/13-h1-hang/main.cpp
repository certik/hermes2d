#include "hermes2d.h"



int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);

  Mesh mesh;
  mesh.load("test.mesh");
  mesh.refine_element(1);
  mesh.refine_element(7, 2);
  mesh.refine_element(2);
  mesh.refine_element(12);
  
  /*for (int i = 0; i < 10; i++)
  {
    Element* e;
    do
    {
      int id = rand() % mesh.get_max_element_id();
      e = mesh.get_element(id);
    }
    while (!e->active);
    
    int ref = 0;
    if (e->is_quad()) ref = rand() % 2 + 1;
    mesh.refine_element(e->id, ref);
  }*/
    
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  H1Space space(&mesh, &shapeset);
  space.set_uniform_order(3);
  space.assign_dofs();
  
  BaseView bv;
  bv.show(&space);
  
  hermes2d_finalize();
  return 0;
}

