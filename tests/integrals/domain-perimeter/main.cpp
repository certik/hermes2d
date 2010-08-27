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
//  The latter may require an additional mesh file to be named adomain.meshX 
//  where X should an integer following the largest already available. The
//  header in the existing mesh files should be helpful in this case.

// If you add a new domain case, put its name here. In addition, create a 
// domain.meshX file and read the existing files to help you create yours.

enum Cases { SQUARE, QUAD, ANNULUS };

//******************************************************************************
// Controls
//
//  The following parameters can be changed:
//

const int P_INIT = 1;        // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 2;  // Number of initial uniform mesh refinements.
const Cases CASE = SQUARE;  // pick the case from the enum above

//******************************************************************************
// Helper functions

//------------------------------------------------------------------------------
// Boundary condition types
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
double BdryLength(MeshFunction* meshfn, int bdryMarker)
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
//      Get spatial edge tangent data
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

} // end of BdryLength()


//******************************************************************************
// Main
//
int main(int argc, char* argv[])
{
 // Load the mesh.
 Mesh mesh;
 H2DReader mloader;

 switch (CASE)
 {
  case SQUARE:
   mloader.load("domain.mesh0", &mesh);
   break;
  case QUAD:
   mloader.load("domain.mesh1", &mesh);
   break;
  case ANNULUS:
   mloader.load("domain.mesh2", &mesh);
   break;
  default: 
   mloader.load("domain.mesh0", &mesh);
 }

 // Perform initial mesh refinements.
 for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();


 // Create H1 space with a default shapeset.
 H1Space space(&mesh, bc_types, NULL, P_INIT);

 Solution sln;
 sln.set_zero(&mesh);

 double perimeter = 0.0;
 double exactPerimeter = 0.0;
 double l1,l2,l3,l4 = 0.0;
 const double pi = 4.0 * atan(1.0);

 switch (CASE)
 {
  case SQUARE:
   // Calculate the length of the four boundaries segments.
   l1 = BdryLength(&sln, 1);
   info("Length of bdry 1 = %g\n", l1);

   l2 = BdryLength(&sln, 2);
   info("Length of bdry 2 = %g\n", l2);

   l3 = BdryLength(&sln, 3);
   info("Length of bdry 3 = %g\n", l3);

   l4 = BdryLength(&sln, 4);
   info("Length of bdry 4 = %g\n", l4);
  
   perimeter = l1+l2+l3+l4;
   info("Computed perimeter = %g\n", perimeter);

   // Set exact value from domain.mesh0 file
   exactPerimeter = 8.0;
   info("Exact perimeter = %g\n", exactPerimeter);

   break;

  case QUAD:
   // Calculate the length of the four boundaries segments.
   l1 = BdryLength(&sln, 1);
   info("Length of bdry 1 = %g\n", l1);

   l2 = BdryLength(&sln, 2);
   info("Length of bdry 2 = %g\n", l2);

   l3 = BdryLength(&sln, 3);
   info("Length of bdry 3 = %g\n", l3);

   l4 = BdryLength(&sln, 4);
   info("Length of bdry 4 = %g\n", l4);
  
   perimeter = l1+l2+l3+l4;
   info("Computed perimeter = %g\n", perimeter);

   // Set exact value from domain.mesh1 file
   exactPerimeter = 3.0 + sqrt(5.0) + sqrt(2.0);
   info("Exact perimeter = %g\n", exactPerimeter);

   break;

  case ANNULUS:
   // Calculate the length of the four boundaries segments.
   l1 = BdryLength(&sln, 1);
   info("Length of bdry 1 = %g\n", l1);

   l2 = BdryLength(&sln, 2);
   info("Length of bdry 2 = %g\n", l2);

   perimeter = l1+l2;
   info("Computed perimeter = %g\n", perimeter);

   // Set exact value from domain.mesh1 file
   exactPerimeter = 2.0 * PI * 1.7;
   info("Exact perimeter = %g\n", exactPerimeter);

   break;

  default:
   cout << " A CASE has not been set " << endl;
   cout << " Bailing out. " << endl; exit(1);
 }
 
 // Display the mesh.
 // (100, 0) is the upper left corner position, 600 x 500 is the window size
 MeshView mview("Mesh", 100, 0, 600, 500);
 mview.show(&mesh);

 // Wait for the view to be closed.
 View::wait();

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

 if (fabs(perimeter - exactPerimeter) < 1e-6) {
   printf("Success!\n");
   return ERROR_SUCCESS;
 }
 else {
   printf("Failure!\n");
   return ERROR_FAILURE;
 }
 

 return 0;

}
