#include "hermes2d.h"

// This example shows how to load a mesh, perform various types
// of "manual"  element refinements, and use keyboard and mouse
// controls.

static char text[] = "\
Click into the image window and:\n\
  press 'm' to show element numbers,\n\
  press 'b' to toggle boundary markers,\n\
  enlarge your window and press 'c' to center the mesh,\n\
  zoom into the mesh using the right mouse button\n\
  and move the mesh around using the left mouse button\n\
    -- in this way you can read the numbers of all elements,\n\
  press 'c' to center the mesh again,\n\
  press 'm' to hide element numbers,\n\
  press 's' to save a screenshot in bmp format\n\
    -- the bmp file can be converted to any standard\n\
       image format using the command \"convert\",\n\
  press 'q' to quit.\n\
  Also see help - click into the image window and press F1.\n";


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // perform some sample initial refinements
  mesh.refine_all_elements();          // Refines all elements.
  mesh.refine_towards_vertex(3, 4);    // Refines mesh towards vertex #3 (4x).
  mesh.refine_towards_boundary(2, 4);  // Refines all elements along boundary 2 (4x).
  mesh.refine_element(86, 0);          // Refines element #86 isotropically.
  mesh.refine_element(112, 0);         // Refines element #112 isotropically.
  mesh.refine_element(84, 2);          // Refines element #84 anisotropically.
  mesh.refine_element(114, 1);         // Refines element #114 anisotropically.

  // display the mesh
  MeshView mview("Hello world!", 100, 100, 500, 500);  // (100, 100) is the upper left corner position
  mview.show(&mesh);                                   // 500 x 500 is the window size

  // practice some keyboard and mouse controls
  printf("%s", text);

  // wait for a view to be closed
  View::wait();
  return 0;
}
