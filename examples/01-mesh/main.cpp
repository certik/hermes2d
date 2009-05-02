#include "hermes2d.h"

// This example shows you how to load a mesh, perform various types
// of initial refinements, and use keyboard and mouse controls. 

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");

  // perform some sample initial refinements
  mesh.refine_element(1);              // refines element #1
  mesh.refine_element(2);              // refines element #2
  mesh.refine_all_elements();          // refines all elements
  mesh.refine_towards_boundary(4, 3);  // refines all elements along boundary 4,
                                       // this is repeated three times
  
  // display the mesh
  MeshView mview("Hello world!", 100, 100, 500, 500);
  mview.show(&mesh);
  
  // practice some keyboard and mouse controls
  printf("Click into the image window and:\n");
  printf("  press 'm' to show element numbers\n");
  printf("    -- note that indices start high since some refinements took place already,\n");
  printf("  resize your window and press 'c' to center the mesh,\n");
  printf("  zoom into the corner using the right mouse button\n"); 
  printf("    -- now you can read the numbers of small elements,\n");
  printf("  move the mesh around using the left mouse button,\n");
  printf("  press 'c' to center the mesh again,\n");
  printf("  press 'm' to hide element numbers,\n");
  printf("  press 's' to save a screenshot,\n");
  printf("  press 'q' to finish.\n");
  View::wait();
  return 0;
}
