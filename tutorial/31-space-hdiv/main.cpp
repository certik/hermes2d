#include "hermes2d.h"
#include "solver_umfpack.h"

// This example shows how to use the Hdiv space and
// visualize finite element basis functions. Note that 
// higher-order basis functions in this space comprise 
// edge functions associated with mesh edges (normal 
// component is zero on the boundary of the element patch
// associated with the edge), and bubble functions 
// associated with elements (normal component is 
// zero on the element boundary).

int P_INIT = 3;

int main(int argc, char* argv[])
{
  if (argc < 2) error("Missing mesh file name parameter.");

  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load(argv[1], &mesh);

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
  int ndof = assign_dofs(&space);

  // visualise the FE basis
  VectorBaseView bview;
  bview.show(&space);

  // wait for all views to be closed
  View::wait();
  return 0;
}

