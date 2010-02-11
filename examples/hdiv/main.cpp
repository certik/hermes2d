#include "hermes2d.h"
#include "solver_umfpack.h"

// This example shows how to use the Hdiv space and
// visualize finite element basis functions

int P_INIT = 3;

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // uniform mesh refinements
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  // initialize the shapeset and the cache
  HdivShapeset shapeset;

  // create the Hdiv space
  HdivSpace space(&mesh, &shapeset);

  // set uniform polynomial degrees
  space.set_uniform_order(P_INIT);

  // enumerate basis functions
  space.assign_dofs();

  // visualise the FE basis
  VectorBaseView bview;
  bview.show(&space);

  // wait for keyboard or mouse input
  View::wait("Waiting for all views to be closed.");
  return 0;
}

