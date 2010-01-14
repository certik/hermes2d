#include "hermes2d.h"

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    printf("please input as this format: <meshfile> \n");
    return ERROR_FAILURE;
  }

  char *file_name = argv[1];

  // load the mesh file
  Mesh mesh;
  ExodusIIReader mloader;
  if (mloader.load(file_name, &mesh)) {
    int ne = mesh.get_num_elements();
    printf("Elements = %d\n", ne);
    for (int eid = 0; eid < ne; eid++) {
	  Element *e = mesh.get_element(eid);
	  printf(" #%d:", e->id);
	  for (int iv = 0; iv < e->nvert; iv++) {
        printf(" %d", e->vn[iv]->id);
	  }
      printf(" | %d\n", e->marker);
	}

    int im = 0;
    printf("Markers\n");
    for (int eid = 0; eid < ne; eid++) {
	  Element *e = mesh.get_element(eid);
	  if (!e->active) continue;

	  int nv = e->nvert;
	  for (int iv = 0; iv < nv; iv++) {
	    Node *nd = e->en[iv];
		if (nd->type == 1 && nd->bnd == 1) { // edge node
		  printf(" %d, %d | %d\n", e->vn[iv]->id, e->vn[(iv + 1) % nv]->id, nd->marker);
		}
	  }
	}

    return ERROR_SUCCESS;
  }
  else {
    printf("failed\n");
    return ERROR_FAILURE;
  }

  printf("Success!\n");
  return ERROR_SUCCESS;
}

