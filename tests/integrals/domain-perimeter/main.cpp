#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#define PI 4.0*atan(1.0)
#include "hermes2d.h"

#include <iostream>
using std::cout;
using std::endl;

//******************************************************************************
//  This is a very basic test for computing the perimeter of a domain. 
//  The user can experiment with existing cases or contribute a new case.
//  The latter may require an additional mesh file to be created and named 
//  domain.meshX where X should be an integer following the largest already 
//  used in the provided set of domain.meshX files in this directory. The
//  header in the existing mesh files should be helpful.

//  To add more domains to test, edit the CMakeLists.txt file. Note: The 
//  boundary markers must be 1, 2, 3, or 4, others will be skipped.

//******************************************************************************
// Controls
//
//  The following parameters can be changed:

const int P_INIT = 1;       // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 2; // Number of initial uniform mesh refinements.

//******************************************************************************
// Helper functions

//------------------------------------------------------------------------------
// Boundary condition types; this does not really affect the testing
// Set your bc_types() according to the CASE
//
BCType bc_types(int marker)
{
  if      (marker == 1) return BC_ESSENTIAL; 
  else if (marker == 2) return BC_NATURAL; 
  else if (marker == 3) return BC_ESSENTIAL; 
  else if (marker == 4) return BC_NATURAL; 
  else                  return BC_NATURAL;
} // end of bc_types()

//------------------------------------------------------------------------------
// Compute marked boundary length 
//
double CalculateBoundaryLength(MeshFunction* meshfn, int bdryMarker)
{
  Quad2D* quad = &g_quad_2d_std;
  meshfn->set_quad_2d(quad);

  Element* e;
  Mesh* mesh = meshfn->get_mesh();

  double length = 0;
  for_all_active_elements(e, mesh) // vfda:  how about something like
  {                                // for_all_active_bdry_elements(e,mesh,bdryMrk)
    for(int edge = 0; edge < e->nvert; ++edge) 
    {
      if ((e->en[edge]->bnd) && (e->en[edge]->marker == bdryMarker)) 
      {
        update_limit_table(e->get_mode());
        RefMap* mu = meshfn->get_refmap();

        meshfn->set_active_element(e); // vfda: should this be here? 
        int eo = quad->get_edge_points(edge);
        meshfn->set_quad_order(eo, H2D_FN_VAL); // vfda: is this needed? 

//      Need quadrature weights...
//      Get quadrature points data for this element order:
//      pt[i][0], pt[i][1] are coordinates; p[i][2] are quadrature weights

        double3* pt  = quad->get_points(eo);

//      Need the magnitude of the spatial tangent vector
//      Get spatial edge tangent data:
//      tan[i][0], tan[i][1] are the x,y components; tan[i][2] is the norm

        double3* tan = mu->get_tangent(edge);

        for (int i = 0; i < quad->get_num_points(eo); i++)
        {
          length += pt[i][2] * tan[i][2];
        }
      }
    }
  }
  
  return 0.5 * length; 

} // end of CalculateBoundaryLength()

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

//******************************************************************************
// Main
//
int main(int argc, char* argv[])
{
  if (argc != 3)
  {
    printf("Missing mesh filename and domain length as command-line parameters.");
    return ERROR_FAILURE;
  }

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load(argv[1], &mesh);

/* Graphics not allowed on test versions; user may uncomment for debugging; vfda
  // Display the mesh.
  // (100, 0) is the upper left corner position, 600 x 500 is the window size
  MeshView mview("Mesh", 100, 0, 600, 500);
  mview.show(&mesh);

  // Wait for the view to be closed.
  View::wait();
*/

  double length = atof(argv[2]);

  // Perform initial mesh refinements.
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create H1 space with a default shapeset.
  H1Space space(&mesh, bc_types, NULL, P_INIT);

  // FIXME: This Solution is artificial here and it should be removed. The 
  // function CalculateBoundaryLength() should only take a pointer to Mesh and 
  // a boundary marker as parameters.
  Solution sln;
  sln.set_zero(&mesh);

  // Calculate the length of the four boundaries segments.
  double l1 = CalculateBoundaryLength(&sln, 1);
  info("Length of boundary 1 = %g\n", l1);

  double l2 = CalculateBoundaryLength(&sln, 2);
  info("Length of boundary 2 = %g\n", l2);

  double l3 = CalculateBoundaryLength(&sln, 3);
  info("Length of boundary 3 = %g\n", l3);

  double l4 = CalculateBoundaryLength(&sln, 4);
  info("Length of boundary 4 = %g\n", l4);
  
  double perimeter = l1 + l2 + l3 + l4;
  info("Computed perimeter = %10.15f\n", perimeter);

  // Set exact value from CMakeLists.txt file
  info("Exact perimeter = %g\n", length);

  if (fabs(perimeter - length) < 1e-6) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
 

  return 0;
}
