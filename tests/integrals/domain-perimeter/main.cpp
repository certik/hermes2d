#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

#include <iostream>
using std::cout;
using std::endl;

//  This is an very basic test for computing the perimeter of a domain. 
//
//  The following parameters can be changed:

const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;       // Number of initial uniform mesh refinements.

// Boundary condition types
BCType bc_types(int marker)
{
  if      (marker == 1) return BC_ESSENTIAL; 
  else if (marker == 2) return BC_NATURAL; 
  else if (marker == 3) return BC_ESSENTIAL; 
  else if (marker == 4) return BC_NATURAL; 
  else                  assert(1);
}

// compute marked boundary length 
double BdryLength(MeshFunction* meshfn, int bdryMarker)
{
  Quad2D* quad = &g_quad_2d_std;
  meshfn->set_quad_2d(quad);

  Element* e;
  Mesh* mesh = meshfn->get_mesh();

  // Show the mesh.
  //MeshView mv("", 0, 0, 400, 400);
  //mv.show(mesh);
  //View::wait();

  int nActiveEle = 0;
  double length = 0;
  for_all_active_elements(e, mesh)         // there should be something like 
  {                                        // for_all_active_bdry_elements(e,mesh)

    ++nActiveEle;
    cout << "ele #               = " << nActiveEle << endl;
    cout << "e->nvert            = " << e->nvert << endl;
    for(int edge = 0; edge < e->nvert; edge++) {
      cout << "e->en[edge]->marker = " << e->en[edge]->marker << endl;
      if ((e->en[edge]->bnd) && (e->en[edge]->marker == bdryMarker)) {
        update_limit_table(e->get_mode());
        RefMap* ru = meshfn->get_refmap();

        meshfn->set_active_element(e);
        int eo = quad->get_edge_points(edge);
        meshfn->set_quad_order(eo, H2D_FN_VAL);


        double3* pt  = quad->get_points(eo);

        for (int i = 0; i < quad->get_num_points(eo); i++) {
          length += pt[i][2];
        }
      }
    }
  }
  
  return 0.5 * length; 

} // end of BdryLength()


int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create H1 space with a default shapeset.
  H1Space space(&mesh, bc_types, NULL, P_INIT);

  Solution sln;
  sln.set_zero(&mesh);

  // Calculate the length of the four boundaries.
  double l1 = BdryLength(&sln, 1);
  info("Length of bdry 1 = %g\n", l1);

  double l2 = BdryLength(&sln, 2);
  info("Length of bdry 2 = %g\n", l2);

  double l3 = BdryLength(&sln, 3);
  info("Length of bdry 3 = %g\n", l3);

  double l4 = BdryLength(&sln, 4);
  info("Length of bdry 4 = %g\n", l4);

  return 0;
}
