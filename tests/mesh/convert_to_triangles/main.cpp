#include <vector>
#include "hermes2d.h"

// This example shows how to convert all elements into triangle elements.

bool do_test(const char* filename) {
  printf("Testing file: %s\n", filename);

  // load the mesh file
  Mesh mesh;
  Element* e;
  mesh.load(filename);

  // Calculate the number of elements after refinement, starting from 0
  int element_num = 0;
  printf("The number of base elements is %d\n", mesh.get_max_element_id());
  for_all_elements(e, &mesh)
  {
    printf("e->id = %d  ", e->id);
    if (e->is_quad())
    {
      printf("type : quadrangle ");
      // one quad element refined into two elements
      element_num += 2;
    }
    else
    {
      printf("type : triangle   ");
      // one triangle element refined into four elements
      element_num += 4;
    }
    if (e->is_curved())
      printf("  curved\n");
    else
      printf("\n");
  }

  // convert the mesh
  mesh.convert_to_triangles();
  if (element_num != mesh.get_max_element_id())
    return false;

  printf("The number of refined elements is %d\n", mesh.get_max_element_id());
  for_all_elements(e, &mesh)
  {
    printf("e->id = %d  ", e->id);
    if (e->is_quad())
    {
      printf("type : quadrangle ");
    }
    else
    {
      printf("type : triangle   ");
    }
    if (e->is_curved())
      printf("  curved\n");
    else
      printf("\n");
  }
  return true;
}

int main(int argc, char* argv[])
{
  printf("test: converting to triangles\n");

  //gather resultss
  const char* test_files[] = {"domain.mesh", "square.mesh", "square_tri.mesh", NULL };
  std::vector<bool> test_results;
  int inx = 0;
  while (test_files[inx] != NULL) {
    test_results.push_back(do_test(test_files[inx]));
    inx++;
  }

  //print results
  printf("\ntest results:\n");
  inx = 0;
  while (test_files[inx] != NULL) {
    if (test_results[inx])
      printf("  test file \"%s\": ok\n", test_files[inx]);
    else
      printf("! test file \"%s\": FAILED\n", test_files[inx]);
    inx++;
  }

  return 0;
}

