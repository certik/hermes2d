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

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load(argv[1], &mesh);

  // Uniform mesh refinements.
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  // Initialize the shapeset and the cache.
  HdivShapeset shapeset;

  // Create an Hdiv space.
  HdivSpace space(&mesh, &shapeset);

  // Set uniform polynomial degrees.
  space.set_uniform_order(P_INIT);

  // Enumerate basis functions.
  int ndof = assign_dofs(&space);

  // Visualise the FE basis.
  VectorBaseView bview;
  bview.show(&space);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

