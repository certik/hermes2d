#include "hermes2d.h"


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");
  mesh.refine_element(1);

  // create a shapeset and an H1 space
  H1Shapeset shapeset;
  H1Space space(&mesh, &shapeset);

  // assign element orders and initialize the space
  space.set_uniform_order(5);
  space.assign_dofs();

  // view the basis functions
  BaseView bview;
  bview.show(&space);

  View::wait();
  return 0;
}

