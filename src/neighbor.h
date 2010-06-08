#ifndef NEIGHBOR_H_
#define NEIGHBOR_H_

#include "common.h"
#include "mesh.h"
#include "quad.h"
#include "solution.h"
#include "forms.h"


/*!\class NeighborSearch neighbor.h "src/neighbor.h"
 * \brief This class serves to find and offer all information from neighbors at a given(active) element and its concrete edge.
 *
 * How works the finding neighbors
 * We will call the active  element as a central element.
 * If we have irregular mesh, there can be three options in relation of central element and neighbor element at common edge.
 * First, the neighbor is same "size" as central, so the edge has active elements on both sides. This option is tested by function get_neighbor().
 * Second, the neighbor is "bigger" then central. Then we have to go "way up" and use method finding_act_elem_up().
 * Third, the neighbor is "smaller", we have more neighbors against the edge. This solves "way down" by method finding_act_elem_down().
 * The choice between way up or way down is made by testing if we can find vertex in the middle of the common edge. If we can
 * then we go way down.

 * Also at every way we fill function values, derivatives and etc. of central and neighbor elements threw method set_fn_values(). Last step is
 * possible change of order of neighbor's function values to correspond function values of central element at same points.

 * We also need transform solution either on neighbor or central element to get points at correct part of the edge. We use method "push_transform"
 * and use only range [0-3]. These types of transformation are common for triangles and quads and choosing right transformation can
 * be derived from local numbers of edges.

 * For numbering and ordering of edges, vertices and sons of an element look into mesh.cpp
 */


class H2D_API NeighborSearch
{
public:

	/*! Common constructor. The least what user has to provide is an active element and mesh.
	* If he wants also function values, he must provide a MeshFunction (solution).
	*  The space is for improvement of the algorithm for choosing the order on the edge. Both, solution and space, have to be defined over the given mesh.
	*/
	NeighborSearch(Element* e, Mesh* mesh, MeshFunction* sln = NULL, Space* space = NULL);

	~NeighborSearch();

	//! Set active edge and compute all needed informations from neighbors.
	// \param[in] edge This is local (element dependent) number of the edge.
	void set_active_edge(int edge);

	//! Number of neighbor elements on edge
	int get_number_of_neighbs();

	//! Return array of transformations of neighbor or central element.
	// \param[in] part_edge Corresponding to the number of neighbor in direction we are enumerate them.
	int* get_transformations(int part_edge);

	//! Return function values of neighbor at integration points.
	// \param[in] part_edge Corresponding to the number of neighbor in direction we are enumerate them.
	scalar* get_fn_values_neighbor(int part_edge);

	//! Return function values of central element at integration points.
	// \param[in] part_edge Corresponding to the number of neighbor in direction we are enumerate them.
	scalar* get_fn_values_central(int part_edge);

	//! Return pointer to function which contains all information from central.
	// \param[in] part_edge Corresponding to the number of neighbor in direction we are enumerate them.
	Func<scalar>* get_values_central(int part_edge);

	//! Return pointer to function which contains all information from neighbor.
	// \param[in] part_edge Corresponding to the number of neighbor in direction we are enumerate them.
	Func<scalar>* get_values_neighbor(int part_edge);


	//! Return number of integration points.
	// \param[in] part_edge Corresponding to the number of neighbor in direction we are enumerate them.
	int get_n_integ_points(int part_edge);

	//! Return pointer to the vector of neighbors id.
	std::vector<int>* get_neighbors();

	//! Return local number of neighbor edge.
	// \param[in] part_edge Corresponding to the number of neighbor in direction we are enumerate them.
	int get_number_neighb_edge(int part_edge);

	//! Return orientation of neighbor edge.
	//! \param[in] part_edge Corresponding to the number of neighbor in direction we are enumerate them.
	//! \return Return 0 if the orientation of the common edge is same on the neighbor element as on the central element.
	//! \return If not than the returned integer is equal to 1.
	int get_orientation_neighb_edge(int part_edge);

	//! Returns a vector where each member contains maximum of orders of central and neighbor element.
	std::vector<int>* get_orders();

	//! Set the order of the integration.
	//! \param[in] order The order for which you want get function values at integration points.
	void set_order_of_integration(int order);

	//! Set MeshFunction (solution), this method is used in case you want to get values of different solutions over the same edge.
	//! Don't have to make another instance of the class.

	void set_solution(MeshFunction* solution);


private:
	const static int max_n_trans = 20;    //!< Number of allowed transformations, see "push_transform" in transform.h.

	int n_neighbors; //!< Number of neighbors.
	Quad2D* quad;
	Mesh* mesh;
	Element* central_el; //!< Central element.
	Element* neighb_el;  //!< Actual neighbor element we are working with,
	MeshFunction* sol;
	Space* space;
	int transformations[max_n_trans][max_n_trans];	//!< Table of transformations for all neighbors.
	int n_trans[max_n_trans];  //!< Number of transformations for every neighbor.
	int active_edge;			     //!< Edge where we are searching for neighbors.
	int neighbor_edge;		   	//!< Edge of the working neighbor corresponding to active_edge.
	scalar* fn_values[max_n_trans]; //!< Function values for central element.
	scalar* fn_values_neighbor[max_n_trans]; //!< Function values for all neighbor elements.
	int np[max_n_trans];						//!< Number of integration points for every neighbor.
	int central_order;  //!< Order of the MeshFunction on central element.
	int neighbor_order; //!< Order of the MeshFunction on the working neighbor element.

	/*! \var
	 * This has two meanings. First it is a flag if the user for setting the order of integration used method set_order_of_integration().
	 * Initial, in the constructor, set equal to -1, else is equal to the order get from the method.
	 */
	int max_of_orders;
	int way_flag; //!< This flag holds which way was used on the active edge. So is equal to one of the members of Trans_flag.



	std::vector<int> neighbors_id; //!<  Vector containing id's of all neighbors. (obsolete)
	std::vector<Element*> neighbors; //!<  Vector containing pointers to  all neighbors.
	std::vector<int> orders; //!< Each member of this vector contains maximum of orders of central and neighbor element.

	//! Vectors of all values (function, derivatives, etc.) for central and neighbor elements.
	std::vector<Func<scalar>*> values_central;
	std::vector<Func<scalar>*> values_neighbor;

	//! Method "way up" for finding neighbor element, from smaller to larger.
	/*!
	 * \param[in] elem The pointer to parent element of the element from previous step.
	 * \param[in] edge_num The active edge.
	 * \param[in] orig_vertex_id Array containing oriented vertices of the active edge.
	 * \param[in] road_vertices Array of vertices which we used in finding the active neighbor.
	 * \param[in] n_road_vertices Number of vertices in array road_vertices.
	 */
	void finding_act_elem_up( Element* elem, int edge_num, int* orig_vertex_id, Node** road_vertices, int n_road_vertices);

	//! Method "Way down" for finding neighbor elements, from larger to smaller.
	/*!
	 * \param[in] vertex The pointer to a vertex which was in the middle of the edge we were working with in previous step.
	 * \param[in] par_vertex_id Array containing id of vertices which are "parents" of the middle vertex (first parameter). They are passed separated because we need to conserve orientation of the vertices.
	 * \param[in] road Array which serves for storing codes of transformations (code is equal to number of son on which we will transform solution).
	 * \param[in] n_road Number of valid members in array road.
	 * \param[in] use_edge Just local name for active edge.
	 *  \param[in] n_vert Number of vertices of the central element.
	 */
	void finding_act_elem_down( Node* vertex, int* par_vertex_id, int* road, int n_road, int use_edge, int n_vert);

	//! Setting the sequence of function values of neighbor in same direction as on central element.
	void set_correct_direction(int index);

	//! Set order of the working edge. Depends on if the space is given.
	int get_max_order();

  /*! Method for getting an order on the working edge. Originally taken from class Space.
   *	\param[in] e Element with which we are working.
   *	\param[in] edge Local edge number of the element e.
   * \return The order of the edge.
   */
	int get_edge_order(Element* e, int edge);

	/*! Helping method for get_edge_order().
   * \param[in] en Edge node whose order we want to find.
   * \return The order of the edge.
   */
	int get_edge_order_internal(Node* en);

  /*! Just reverse values in vector.
   * \param[in] vector The pointer to array in which we will reverse values.
   * \param[in] n Length of the vector.
   */
	void reverse_vector(scalar* vector, int n);

  //! Structure containing all needed information about neighbor's edge. Initial values of both members are invalid.

	struct NeighborEdgeInfo{
		NeighborEdgeInfo(){
			local_num_of_edge = -1;
			orientation = -1;
		}

		//! Local number of the edge on neighbor element.
		int local_num_of_edge;

   /*! Relative orientation of the neighbor edge. If equal to 0 then the neighbor edge has same orientation as the active edge,
    * otherwise is equal to 1.
    */
		int orientation;
	};

	/*! Find the orientation of the neighbor edge in relation with the active edge. All input paramaters depend on in which way we
   * we are. For all ways parent1 and parent2 are ids of vertices oriented in same direction as are vertices of central element.
   * Parameter part_of_edge is 0 or 1 and actualy is used only in way down. In other ways is set to 0.
   * For way up parameters parent1 and parent2 are vertices which define the edge of inactive parent of central element. The edge has on other
   * side active neighbor.
   * For way down parent1 and parent2 are vertices of the edge, which is part of active edge. Both vertices belong to direct parent element
   * of active neighbor element. This means that one of them for sure belongs to active neighbor. The second vertex serves for finding the other(middle) vertex
   * which define neighbor edge corresponding to active edge.
   * Parameter part_of_edge is 0 if the neighbor edge of active neighbor has vertices parent1 and middle vertex and equal to 1 if the vertices which define
   * the neighbor edge are middle vertex and parent2.
   * \param[out] edge_info The relative orientation is copied into the struct NeighborEdgeInfo.
   */
	void direction_neighbor_edge(int parent1, int parent2, int part_of_edge, NeighborEdgeInfo* edge_info);

	// cleaning before usage of given edge
	void clean_all();

	/*! This serves for distinguish which way was used for finding neighbors and according the way how are obtain values in method
   * set_fn_values().
   */
	enum Trans_flag{
		H2D_NO_TRANSF = 0,  //!< Don't use any transformation, the edge has on both sides active element.
		H2D_WAY_DOWN = 1, 	//!< Transformation of central element, against the edge neighbor has some sons.
		H2D_WAY_UP = 2			//!< Transformation of neighbor element, central element is son.
	};

	/*! Fill function values of central and neighbor element.
   * \param[in] flag Flag by which we decide on what element will be applicated transformations.
   */
	void set_fn_values(Trans_flag flag);

	//! Vector containing all neighbor edges information corresponding to active edge.
	std::vector<NeighborEdgeInfo> neighbor_edges;

	//! This method serves for fill all values of central and all neighbors elements.
	void compute_fn_values();

};


#endif /* NEIGHBOR_H_ */
