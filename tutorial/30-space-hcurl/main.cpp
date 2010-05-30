#include "hermes2d.h"
#include "solver_umfpack.h"

// This example shows how to use the Hcurl space and
// visualize finite element basis functions. Note that 
// higher-order basis functions in this space comprise 
// edge functions associated with mesh edges (tangential 
// component is zero on the boundary of the element patch
// associated with the edge), and bubble functions 
// associated with elements (tangential component is 
// zero on the element boundary).
//
// The following parameters can be changed:

int INIT_REF_NUM = 2;      // Initial uniform mesh refinement.
int P_INIT = 3;            // Polynomial degree of mesh elements.

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinement.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize the shapeset and the cache.
  HcurlShapeset shapeset;

  // Create an Hcurl space.
  HcurlSpace space(&mesh, &shapeset);

  // Set uniform polynomial degrees.
  space.set_uniform_order(P_INIT);

  // Enumerate basis functions.
  int ndof = assign_dofs(&space);

  // Visualize FE basis.
  VectorBaseView bview;
  bview.show(&space);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

