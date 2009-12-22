// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __HERMES2D_FUNCTION_H
#define __HERMES2D_FUNCTION_H

#include "common.h"
#include "transform.h"
#include "quad_all.h"


// Precalculation masks
enum
{
  FN_VAL_0 = 0x0001, FN_VAL_1 = 0x0040, // Function values
  FN_DX_0  = 0x0002, FN_DX_1  = 0x0080, // First derivative
  FN_DY_0  = 0x0004, FN_DY_1  = 0x0100, // First derivative
  FN_DXX_0 = 0x0008, FN_DXX_1 = 0x0200, // Second derivative
  FN_DYY_0 = 0x0010, FN_DYY_1 = 0x0400, // Second derivative
  FN_DXY_0 = 0x0020, FN_DXY_1 = 0x0800  // Second mixed derivative
};

// Both components are usually requested together...
const int FN_VAL = FN_VAL_0 | FN_VAL_1;
const int FN_DX  = FN_DX_0  | FN_DX_1;
const int FN_DY  = FN_DY_0  | FN_DY_1;
const int FN_DXX = FN_DXX_0 | FN_DXX_1;
const int FN_DYY = FN_DYY_0 | FN_DYY_1;
const int FN_DXY = FN_DXY_0 | FN_DXY_1;

// The most common case is FN_DEFAULT, ie. values and the gradient.
const int FN_DEFAULT = FN_VAL | FN_DX | FN_DY;            // default precalculation mask
const int FN_ALL = FN_DEFAULT | FN_DXX | FN_DYY | FN_DXY; // precalculate everything

const int FN_COMPONENT_0 = FN_VAL_0 | FN_DX_0 | FN_DY_0 | FN_DXX_0 | FN_DYY_0 | FN_DXY_0;
const int FN_COMPONENT_1 = FN_VAL_1 | FN_DX_1 | FN_DY_1 | FN_DXX_1 | FN_DYY_1 | FN_DXY_1;

// Plenty of checking stuff for the debug version
#ifndef NDEBUG
  #define check_params \
    if (component < 0 || component > num_components) \
      error("Invalid component. You are probably using scalar-valued shapeset for an Hcurl problem."); \
    if (cur_node == NULL) \
      error("Invalid node. Did you call set_quad_order()?");
  #define check_table(n, msg) \
    if (cur_node->values[component][n] == NULL) \
      error(msg " not precalculated for component %d. Did you call set_quad_order() with correct mask?", component)
#else
  #define check_params
  #define check_table(n, msg)
#endif


/// \brief Represents an arbitrary function defined on an element.
///
/// The Function class is an abstraction of a function defined in integration points on an
/// element. You first specify what quadrature tables you want to use (set_quad_2d()) and select
/// an element (Transformable::set_active_element()). Then you select concrete integration points
/// (set_quad_order()) and obtain the function values by calling one of the functions get_fn_values(),
/// get_dx_values(), etc.
///
/// This class is a template for RealFunction and ScalarFunction, depending of which type the
/// function values are. For example, shape functions are always real (see PrecalcShapeset), while
/// the solution can be complex (see Solution).
///
/// The design goal for this class is to define a single common interface for functions used as
/// integrands in the weak formulation. It should not matter whether you are integrating a shape
/// function or, for example, a previous solution of the PDE in time-dependent problems.
/// Ideally, you should also be able to apply the bilinear form not only to shape functions
/// during assembling, but also to the solution when calculating energy norms etc. The last
/// feature is unfortunately limited to real code, because a PDE solution can be complex (hence
/// Solution inherits from ScalarFunction), but shape functions are real and for efficiency
/// the bilinear form only takes RealFunction arguments.
///
/// Since this class inherits from Transformable, you can obtain function values in integration
/// points transformed to sub-areas of the current element (see push_transform(), pop_transform()).
///
template<typename TYPE>
class Function : public Transformable
{
public:

  /// Default constructor.
  Function();
  /// Destructor.
  virtual ~Function();


  /// \brief Returns the polynomial degree of the function being represented by the class.
  int get_fn_order() const { return order; }

  /// \brief Returns the number of components of the function being represented by the class.
  int get_num_components() const { return num_components; }


  /// Activates an integration rule of the specified order. Subsequent calls to
  /// get_values(), get_dx_values() etc. will be returning function values at these points.
  /// \param order [in] Integration rule order.
  /// \param mask [in] A combination of one or more of the constants FN_VAL, FN_DX, FN_DY,
  ///   FN_DXX, FN_DYY, FN_DXY specifying the values which should be precalculated. The default is
  ///   FN_VAL | FN_DX | FN_DY. You can also use FN_ALL to precalculate everything.
  void set_quad_order(int order, int mask = FN_DEFAULT)
  {
    pp_cur_node = (void**) JudyLIns(nodes, order, NULL);
    // if you get SIGSEGV here, you maybe forgot to include the function in the list
    // of external functions in WeakForm::add_biform()...
    cur_node = (Node*) *pp_cur_node;
    // another reason may be a bug in Judy array usage in Hermes2D (this needs to be fixed)
    // -- as a workaround, you may for the time being use more than one pss for problems
    // where the basis and test functions can be on different meshes (i.e., multi-mesh).
    if (cur_node == NULL || (cur_node->mask & mask) != mask) precalculate(order, mask);
  }

  /// \brief Returns function values.
  /// \param component [in] The component of the function (0 or 1).
  /// \return The values of the function at all points of the current integration rule.
  TYPE* get_fn_values(int component = 0)
  {
    check_params; check_table(0, "Function values");
    return cur_node->values[component][0];
  }

  /// \brief Returns the x partial derivative.
  /// \param component [in] The component of the function (0 or 1).
  /// \return The x partial derivative of the function at all points of the current integration rule.
  TYPE* get_dx_values(int component = 0)
  {
    check_params; check_table(1, "DX values");
    return cur_node->values[component][1];
  }

  /// \brief Returns the y partial derivative.
  /// \param component [in] The component of the function (0 or 1).
  /// \return The y partial derivative of the function at all points of the current integration rule.
  TYPE* get_dy_values(int component = 0)
  {
    check_params; check_table(2, "DY values");
    return cur_node->values[component][2];
  }

  /// \brief Returns both x and y partial derivatives.
  /// This function provides the both often-used dx and dy values in one call.
  /// \param dx [out] Variable which receives the pointer to the first partial derivatives by x
  /// \param dy [out] Variable which receives the pointer to the first partial derivatives by y
  /// \param component [in] The component of the function (0 or 1).
  void get_dx_dy_values(TYPE*& dx, TYPE*& dy, int component = 0)
  {
    check_params; check_table(1, "DX values"); check_table(2, "DY values");
    dx = cur_node->values[component][1];
    dy = cur_node->values[component][2];
  }

  /// \brief Returns the second x partial derivative.
  /// \param component [in] The component of the function (0 or 1).
  /// \return The x second partial derivative of the function at all points of the current integration rule.
  TYPE* get_dxx_values(int component = 0)
  {
    check_params; check_table(3, "DXX values");
    return cur_node->values[component][3];
  }

  /// \brief Returns the second y partial derivative.
  /// \param component [in] The component of the function (0 or 1).
  /// \return The y second partial derivative of the function at all points of the current integration rule.
  TYPE* get_dyy_values(int component = 0)
  {
    check_params; check_table(4, "DYY values");
    return cur_node->values[component][4];
  }

  /// \brief Returns the second mixed derivative.
  /// \param component [in] The component of the function (0 or 1).
  /// \return The second mixed derivative of the function at all points of the current integration rule.
  TYPE* get_dxy_values(int component = 0)
  {
    check_params; check_table(5, "DXY values");
    return cur_node->values[component][5];
  }

  /// For internal use.
  TYPE* get_values(int a, int b)
  {
    return cur_node->values[a][b];
  }


  /// \brief Selects the quadrature points in which the function will be evaluated.
  /// \details It is possible to switch back and forth between different quadrature
  /// points: no precalculated values are freed. The standard quadrature is
  /// always selected by default already.
  /// \param quad_2d [in] The quadrature points.
  virtual void set_quad_2d(Quad2D* quad_2d);

  /// \brief Returns the current quadrature points.
  Quad2D* get_quad_2d() const { return quads[cur_quad]; }


  /// See Transformable::push_transform()
  virtual void push_transform(int son);

  /// See Transformable::pop_transform()
  virtual void pop_transform();


  /// \brief Frees all precalculated tables.
  virtual void free() = 0;


protected:

  /// precalculates the current function at the current integration points.
  virtual void precalculate(int order, int mask) = 0;

  int order;          ///< current function polynomial order
  int num_components; ///< number of vector components

  struct Node
  {
    int mask;           ///< a combination of FN_XXX: specifies which tables are present
    int size;           ///< size in bytes of this struct (for maintaining total_mem)
    TYPE* values[2][6]; ///< pointers to 'data'
    TYPE data[1];       ///< value tables. The length may vary.
  private: //operation that are not allowed due to the variable length of the Node structure
    Node(const Node& org) {}; ///< Copy constructor is disabled.
    Node& operator=(const Node& other) { return *this; }; ///< Assignment is not allowed.
  };

  void** sub_tables;  ///< pointer to the current secondary Judy array
  void** nodes;       ///< pointer to the current tertiary Judy array FIXME
  void** pp_cur_node;
  void*  overflow_nodes;
  Node*  cur_node;

  void update_nodes_ptr()
  {
    if (sub_idx > max_idx)
      handle_overflow_idx();
    else { 
      debug_assert((sub_idx >> (sizeof(Word_t) * 8)) == 0, "E index is larger than JudyLins can contain (Function::update_nodes_ptr)");
      nodes = (void**) JudyLIns(sub_tables, (Word_t)sub_idx, NULL);
    }
  }

  /// For internal use only.
  void force_transform(uint64_t sub_idx, Trf* ctm)
  {
    this->sub_idx = sub_idx;
    this->ctm = ctm;
    update_nodes_ptr();
  }

  Quad2D* quads[4]; ///< list of available quadratures
  int cur_quad;     ///< active quadrature (index into 'quads')
  int total_mem;    ///< total memory in bytes used by the tables
  int max_mem;      ///< peak memory usage

  Node* new_node(int mask, int num_points); ///< allocates a new Node structure
  void  free_nodes(void** nodes);
  void  free_sub_tables(void** sub);
  void  handle_overflow_idx();

  void replace_cur_node(Node* node)
  {
    if (cur_node != NULL) { total_mem -= cur_node->size; ::free(cur_node); }
    *pp_cur_node = node;
    cur_node = node;
  }

  void check_order(Quad2D* quad, int order)
  {
    if (order < 0 || order >= quad->get_num_tables())
      error("Order out of range (%d, %d).", order, quad->get_num_tables());
  }

  static int idx2mask[6][2];  ///< index to mask table

};


/// Represents a real function on an element.
typedef Function<double> RealFunction;

/// Represents a scalar function on an element.
typedef Function<scalar> ScalarFunction;


#undef check_params
#undef check_table



//// implementation of non-inline template members /////////////////////////////////////////////////


template<typename TYPE>
Function<TYPE>::Function()
              : Transformable()
{
  order = 0;
  max_mem = total_mem = 0;

  nodes = NULL;
  cur_node = NULL;
  sub_tables = NULL;
  overflow_nodes = NULL;

  memset(quads, 0, sizeof(quads));
}


template<typename TYPE>
Function<TYPE>::~Function()
{
  if (overflow_nodes != NULL)
    free_nodes(&overflow_nodes);
}


template<typename TYPE>
void Function<TYPE>::set_quad_2d(Quad2D* quad_2d)
{
  int i;

  // check to see if we already have the quadrature
  for (i = 0; i < 4; i++)
    if (quads[i] == quad_2d) {
      cur_quad = i;
      return;
    }

  // if not, add the quadrature to a free slot
  for (i = 0; i < 4; i++)
    if (quads[i] == NULL) {
      quads[i] = quad_2d;
      cur_quad = i;
      return;
    }

  error("too many quadratures.");
}


template<typename TYPE>
void Function<TYPE>::push_transform(int son)
{
  Transformable::push_transform(son);
  if (sub_tables) update_nodes_ptr(); // fixme
}


template<typename TYPE>
void Function<TYPE>::pop_transform()
{
  Transformable::pop_transform();
  if (sub_tables) update_nodes_ptr(); // fixme
}


template<typename TYPE>
int Function<TYPE>::idx2mask[6][2] =
{
  { FN_VAL_0, FN_VAL_1 }, { FN_DX_0,  FN_DX_1  }, { FN_DY_0,  FN_DY_1  },
  { FN_DXX_0, FN_DXX_1 }, { FN_DYY_0, FN_DYY_1 }, { FN_DXY_0, FN_DXY_1 }
};


template<typename TYPE>
typename Function<TYPE>::Node* Function<TYPE>::new_node(int mask, int num_points)
{
  // get the number of tables
  int nt = 0, m = mask;
  if (num_components < 2) m &= FN_VAL_0 | FN_DX_0 | FN_DY_0 | FN_DXX_0 | FN_DYY_0 | FN_DXY_0;
  while (m) { nt += m & 1; m >>= 1; }

  // allocate a node including its data part, init table pointers
  int size = offsetof(Node, data) + sizeof(TYPE) * num_points * nt; //Due to impl. reasons, the structure Node has non-zero length of data even though they can be zero.
  Node* node = (Node*) malloc(size);
  node->mask = mask;
  node->size = size;
  memset(node->values, 0, sizeof(node->values));
  TYPE* data = node->data;
  for (int j = 0; j < num_components; j++) {
    for (int i = 0; i < 6; i++)
      if (mask & idx2mask[i][j]) {
        node->values[j][i] = data;
        data += num_points;
      }
  }
  // todo: maybe put here copying of the old node

  total_mem += size;
  if (max_mem < total_mem) max_mem = total_mem;
  return node;
}


template<typename TYPE>
void Function<TYPE>::free_nodes(void** nodes)
{
  // free all nodes stored in the tertiary Judy array
  unsigned long order = 0;
  void** pp = (void**) JudyLFirst(*nodes, &order, NULL);
  while (pp != NULL)
  {
    // free the concrete Node structure
    total_mem -= ((Node*) *pp)->size;
    ::free(*pp);
    pp = JudyLNext(*nodes, &order, NULL);
  }
  JudyLFreeArray(nodes, NULL);
}


template<typename TYPE>
void Function<TYPE>::free_sub_tables(void** sub)
{
  // iterate through the specified secondary (sub_idx) Judy array
  unsigned long idx = 0;
  void** nodes = (void**) JudyLFirst(*sub, &idx, NULL);
  while (nodes != NULL)
  {
    free_nodes(nodes);
    nodes = JudyLNext(*sub, &idx, NULL);
  }
  JudyLFreeArray(sub, NULL);
}


template<typename TYPE>
void Function<TYPE>::handle_overflow_idx()
{
  if (overflow_nodes != NULL) free_nodes(&overflow_nodes);
  overflow_nodes = NULL;
  nodes = &overflow_nodes;
}




#endif
