#include "hermes2d.h"

// This example demonstrates how to create a space over a mesh
// and how to visualize basis functions. After you run this
// example, try to change P_INIT and/or the mesh and try again!

int P_INIT = 5;

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");
  mesh.refine_element(1);    // you can perform arbitrary other initial
                             // refinements here, see example 01

  // create a shapeset and an H1 space
  H1Shapeset shapeset;
  H1Space space(&mesh, &shapeset);

  // assign element orders and initialize the space
  space.set_uniform_order(P_INIT);  // set uniform polynomial order to all mesh elements
  space.assign_dofs();         // enumerate basis functions
                               // this needs to be done whenever the space changes

  // view the basis functions
  BaseView bview;
  bview.show(&space);
  printf("Click into the image window and:\n");
  printf("  press 'f' to make the color scale finer/coarser,\n");
  printf("  press '3' to see 3D plot of the basis functions,\n");
  printf("  use all three mouse buttons to rotate/move/enlarge the graphs,\n");
  printf("  use the right/left arrows to browse through basis functions,\n");
  printf("  press 'l' to see linearizer output\n");
  printf("    -- basis functions are plotted using piecewise-linear approximation (welcome to OpenGL),\n");
  printf("  use the left mouse button to drag the scale to another corner,\n");
  printf("  press 'p' to switch to greyscales and back,\n");
  printf("  press 'm' to hide the mesh,\n");
  printf("  press 's' to save screenshots,\n");
  printf("  press 'q' to finish.\n");

  // wait for keyboard or mouse input
  View::wait();
  return 0;
}

