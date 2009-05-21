#include "hermes2d.h"
#include "solver_umfpack.h"


int P_INIT = 3;

int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("square.mesh");
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  HdivShapeset shapeset;
  HdivSpace space(&mesh, &shapeset);

  space.set_uniform_order(P_INIT);
  space.assign_dofs();

  VectorBaseView bview;
  bview.show(&space);

  View::wait();
  return 0;
}

