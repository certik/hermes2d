#include "hermes2d.h"

// This example demonstrates how to create a space over a mesh
// and how to visualize the finite element basis functions
// using the BaseView class. The variable P_INIT defines the
// initial degree of all mesh elements. Currently, it is set
// to 1. After visualizing the basis functions using the hints
// printed in the terminal window, change the value of P_INIT
// to 2, 3, etc. to also see higher-order shape functions.
//
// You can use this example to visualize all shape functions
// on the reference square and reference triangle domains,
// just load the corresponding mesh at the beginning of the
// main.cpp file.

int P_INIT = 1;

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
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);            // original L-shape domain
  //mloader.load("domain_quad.mesh", &mesh);     // reference square
  //mloader.load("domain_tri.mesh", &mesh);      // reference triangle

  // sample element refinement, to see more basis functions
  //mesh.refine_all_elements();

  // create a shapeset and an H1 space
  H1Shapeset shapeset;
  H1Space space(&mesh, &shapeset);

  // assign element orders and initialize the space
  space.set_uniform_order(P_INIT);  // Set uniform polynomial order
                                    // P_INIT to all mesh elements.
  space.assign_dofs();              // Enumerate basis functions.

  // view the basis functions
  BaseView bview;
  bview.show(&space);

  // practice some keyboard and mouse controls
  printf("%s", text);

  // wait for a view to be closed
  View::wait("Waiting for a view to be closed.");
  return 0;
}

