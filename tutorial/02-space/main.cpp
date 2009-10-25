#include "hermes2d.h"

// This example demonstrates how to create a space over a mesh
// and how to visualize finite element basis functions. Below,
// the variable P_INIT defines the initial degree of all mesh
// elements. Currently, it is set to 3. After visualizing the
// basis functions using the printed hints , try to change the
// value of P_INIT, type "make", and run the executable again.

int P_INIT = 3;

static char text[] = "\
Click into the image window and:\n\
  press 'f' to make the color scale finer/coarser,\n\
  press '3' to see 3D plot of the basis functions,\n\
  use all three mouse buttons to rotate/move/enlarge the graphs,\n\
  use the right/left arrows to browse through basis functions,\n\
  press 'l' to see linearizer output\n\
      -- higher-order polynomials are plotted using \n\
         adaptive piecewise-linear approximation,\n\
  use the left mouse button to drag the scale to another corner,\n\
  press 'p' to switch to greyscales and back,\n\
  press 'm' to hide the mesh,\n\
  press 's' to save screenshots,\n\
  press 'q' to quit.\n\
  Also see help - click into the image window and press F1.\n";


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");
  mesh.refine_element(1);           // Sample "manual" element refinement.

  // create a shapeset and an H1 space
  H1Shapeset shapeset;
  H1Space space(&mesh, &shapeset);

  // assign element orders and initialize the space
  space.set_uniform_order(P_INIT);  // Set uniform polynomial order P_INIT to all mesh elements.
  space.assign_dofs();              // Enumerate basis functions.

  // view the basis functions
  BaseView bview;
  bview.show(&space);

  // practice some keyboard and mouse controls
  printf("%s", text);

  // wait for keyboard or mouse input
  printf("Waiting for keyboard or mouse input.\n");
  View::wait();
  return 0;
}

