#include "hermes2d.h"


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");
  mesh.refine_element(1);
  mesh.refine_element(2);
  
  // display the mesh
  MeshView mview("Hello world!", 100, 100, 500, 500);
  mview.show(&mesh);
  
  // wait until the user closes the MeshView window
  View::wait();
  return 0;
}
