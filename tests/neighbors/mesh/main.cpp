#include "hermes2d.h"

// This test tests if class neighbor finds for every active element all its neighbors

#define H2D_ERROR_SUCCESS   0
#define H2D_ERROR_FAILURE  -1

int main(int argc, char* argv[])
{
	// load mesh file
	Mesh mesh;
	H2DReader mloader;
	mloader.load("domain.mesh", &mesh);

	// perform refinement to get complex mesh
	mesh.refine_element(0, 1);
	mesh.refine_element(2);
	mesh.refine_element(3, 1);
	mesh.refine_element(7, 2);
	mesh.refine_element(13, 2);
	mesh.refine_element(8);
	mesh.refine_element(9);
	mesh.refine_element(15);
	mesh.refine_element(17);
	mesh.refine_element(24);
	mesh.refine_element(28);

  // display the mesh
	// MeshView mview("neighbors_test", 100, 100, 500, 500);
	// mview.show(&mesh);
  // wait for keyboard or mouse input
  // View::wait("Waiting for keyboard or mouse input.");

	Element* e = NULL;
	NeighborSearch* neighb = NULL;
	std::vector<int>* neighbors_id;
	std::map<int, std::vector<int> > all_neighbors;
	std::map<int, std::vector<int> >::iterator it; 

	// for every element we save all its neighbors into "all_neighhbors".
	// For every neighbor, his id is compared with already inserted elements into all_neighbors. If is found then
	// in his own vector of his neighbors id the id of active element is searched. Failure of the test is if
	// is not found.

	int local_neighb;
	int test;

	for_all_active_elements(e, &mesh){
		test = 0;
		neighb = new NeighborSearch(e, &mesh);
		for(int i = 0; i < e->nvert; i++){
			if(e->en[i]->bnd == 0)
				neighb->set_active_edge(i);
		}
		neighbors_id = neighb->get_neighbors();
		all_neighbors[e->id] = *neighbors_id;
		delete neighb;
		for(int i = 0; i < all_neighbors[e->id].size(); i++){
			local_neighb = all_neighbors[e->id][i];
			it = all_neighbors.find(local_neighb);
			if(it != all_neighbors.end()){
				for(int j = 0; j < it->second.size(); j++){
					if(it->second[j] == e->id)
						test = 1;
				}
				if(test == 0){
				  printf("Failure!\n");
			    return H2D_ERROR_FAILURE;
					}
			}
		}
	}

  // if you return this, the test will succeed:
  printf("Success! \n");
	return H2D_ERROR_SUCCESS;
}
