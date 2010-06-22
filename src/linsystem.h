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

#ifndef __H2D_LINSYSTEM_H
#define __H2D_LINSYSTEM_H

#include "common.h"
#include "tuple.h"
#include "matrix.h"
#include "matrix_old.h"
#include "forms.h"
#include "weakform.h"
#include <map>

class Space;
class PrecalcShapeset;
class WeakForm;
class CommonSolver;

// Default H2D projection norm in H1 norm.
extern int H2D_DEFAULT_PROJ_NORM;

/// Instantiated template. It is used to create a clean Windows DLL interface.
H2D_API_USED_TEMPLATE(Tuple<int>);
H2D_API_USED_TEMPLATE(Tuple<Space*>);
H2D_API_USED_TEMPLATE(Tuple<MeshFunction*>);
H2D_API_USED_TEMPLATE(Tuple<Solution*>);
H2D_API_USED_TEMPLATE(Tuple<PrecalcShapeset*>);

/// For projection, the user may provide bi/linear forms for each solution component stored
/// in tuples of following types
typedef Tuple< std::pair<WeakForm::biform_val_t, WeakForm::biform_ord_t> > biforms_tuple_t;
typedef Tuple< std::pair<WeakForm::liform_val_t, WeakForm::liform_ord_t> > liforms_tuple_t;

///
///
///
///
///
class H2D_API LinSystem
{
public:

  LinSystem();
  LinSystem(WeakForm* wf_, CommonSolver* solver_);
  LinSystem(WeakForm* wf_);                  // solver will be set to NULL and default solver will be used
  LinSystem(WeakForm* wf_, CommonSolver* solver_, Space* s_);
  LinSystem(WeakForm* wf_, Space* s_);       // solver will be set to NULL and default solver will be used
  LinSystem(WeakForm* wf_, CommonSolver* solver_, Tuple<Space*> spaces_);
  LinSystem(WeakForm* wf_, Tuple<Space*> spaces_);      // solver will be set to NULL and default solver will be used
  LinSystem(WeakForm* wf_, CommonSolver* solver_, Space* space1_, Space* space2_);

  virtual ~LinSystem();

  void init_lin(WeakForm* wf, CommonSolver* solver);
  void init_spaces(Tuple<Space*> spaces);
  void init_space(Space* s);         // single equation case
  void set_spaces(Tuple<Space*> spaces);
  void set_pss(Tuple<PrecalcShapeset*> pss);
  void set_pss(PrecalcShapeset* p);  // single equation case
  void copy(LinSystem* sys);
  Space* get_space(int n) {
      if (n < 0 || n >= this->wf->neq) error("Bad index of space.");
      return this->spaces[n];
  }
  Mesh* get_mesh(int n) {
      if (n < 0 || n >= this->wf->neq) error("Bad index of mesh.");
      return this->spaces[n]->mesh;
  }
  PrecalcShapeset* get_pss(int n) {
      if (n < 0 || n >= this->wf->neq) error("Bad index of precalc shapeset.");
      return this->pss[n];
  }

  /// Helps to determine if linear or nonlinear class instance is used
  /// similar to Java instance of functionality
  virtual bool is_linear() { return true; }

  /// Assembles the stiffness matrix and load vector. Vectors Vec, Dir and
  /// RHS must be allocated when assemble() is called.
  virtual void assemble(bool rhsonly = false);
  void assemble_rhs_only() { assemble(true); }

  /// Solves the matrix problem and propagates the resulting coefficient vector into
  /// one or more Solutions. The solution class does not contain the original solution
  /// vector. Instead it contains a new coefficient vector that corresponds to a monomial
  /// In this way, Solution does not require a copy of Space. In other words, Solution
  /// contains the last copy of the vector Vec even after this vector is freed in consequent
  /// computation. This is used in algorithms that require previous solutions, such as
  /// the Newton's method, time stepping, etc.
  bool solve(Tuple<Solution*> sln);
  bool solve(Solution* sln); // single equation case
  bool solve(Solution* sln1, Solution* sln2); // two equations case
  bool solve(Solution* sln1, Solution* sln2, Solution* sln3); // three equations case

  /// Frees the stiffness matrix.
  virtual void free_matrix();

  /// Frees the stiffness matrix, coefficient vectors, and matrix solver data.
  virtual void free();

  /// Saves the stiffness matrix in various formats.
  void save_matrix_matlab(const char* filename, const char* varname = "A");
  void save_rhs_matlab(const char* filename, const char* varname = "b");
  void save_matrix_bin(const char* filename);
  void save_rhs_bin(const char* filename);

  void enable_dir_contrib(bool enable = true) {  want_dir_contrib = enable;  }

  scalar* get_solution_vector() { return Vec; }
  int get_num_dofs();
  int get_num_dofs(int i) {
    if (this->spaces[i] == NULL) error("spaces[%d] is NULL in LinSystem::get_num_dofs().", i);
    return this->spaces[i]->get_num_dofs();
  }
  int get_num_spaces() { return this->wf->neq; };
  int get_matrix_size();
  void get_matrix(int*& Ap, int*& Ai, scalar*& Ax, int& size);
  void get_rhs(scalar*& RHS, int& size) { RHS = this->RHS; size=this->get_num_dofs(); }
  void get_solution_vector(scalar*& sln_vector, int& sln_vector_len)
       { sln_vector = Vec; sln_vector_len = this->get_num_dofs(); }

  /// Returns a copy of the solution vector.
  void get_solution_vector(std::vector<scalar>& sln_vector_out);

  /// Allocate vectors Vec, RHS and Dir of length this->ndof. All vectors
  /// must be NULL at input, to make sure that the user is not losing
  /// information stored in these vectors.
  void alloc_and_zero_vectors();

  /// Reallocates vectors Vec, RHS and Dir according to a new this->ndof.
  void realloc_and_zero_vectors();

  /// Frees vectors Vec, RHS and Dir and sets them to NULL. This should
  /// be used very carefully since the vector Vec stores the actual solution
  /// coefficients. Typicaly this needs to be done after the space changes
  /// and thus the vector Vec loses its meaning. Before freeing it, however,
  /// the information contained in it should be recycled, for example via
  /// a projection onto the new space.
  void free_vectors();

  /*
  /// For debug purposes.
  void print_vector();
  */

  /// Assigning DOF = enumerating basis functions in the FE spaces.
  int assign_dofs();  // all spaces

  /// Global orthogonal projection of multiple solution components. For each of
  /// them a different proj_norm can be used. This defines the entire coefficient
  /// vector Vec. Calls assign_dofs() at the beginning.
  void project_global(Tuple<MeshFunction*> source, Tuple<Solution*> target, Tuple<int> proj_norms = Tuple<int>());

  /// The same as above, but the user may specify the forms that are used in the projection
  /// (useful e.g. when working in curvilinear coordinate systems).
  void project_global(Tuple<MeshFunction*> source, Tuple<Solution*> target, biforms_tuple_t proj_biforms, liforms_tuple_t proj_liforms);

  /// Global orthogonal projection of one MeshFunction.
  void project_global(MeshFunction* source, Solution* target, int proj_norm = H2D_DEFAULT_PROJ_NORM)
  {
    if (this->wf->neq != 1)
      error("Number of projected functions must be one if there is only one equation, in LinSystem::project_global().");
    this->project_global(Tuple<MeshFunction*>(source), Tuple<Solution*>(target), Tuple<int>(proj_norm));
  };

  /// Global orthogonal projection of one MeshFunction -- user specified projection bi/linear forms.
  void project_global(MeshFunction* source, Solution* target,
                  std::pair<WeakForm::biform_val_t, WeakForm::biform_ord_t> proj_biform,
                  std::pair<WeakForm::liform_val_t, WeakForm::liform_ord_t> proj_liform)
  {
    if (this->wf->neq != 1)
      error("Number of projected functions must be one if there is only one equation, in LinSystem::project_global().");
    this->project_global(Tuple<MeshFunction*>(source), Tuple<Solution*>(target), biforms_tuple_t(proj_biform), liforms_tuple_t(proj_liform));
  };


  /// Global orthogonal projection of one scalar ExactFunction.
  void project_global(ExactFunction source, Solution* target, int proj_norm = H2D_DEFAULT_PROJ_NORM)
  {
    if (this->wf->neq != 1)
      error("Number of projected functions must be one if there is only one equation, in LinSystem::project_global().");
    if (proj_norm != 0 && proj_norm != 1) error("Wrong norm used in orthogonal projection (scalar case).");
    Mesh *mesh = this->get_mesh(0);
    if (mesh == NULL) error("Mesh is NULL in project_global().");
    Solution sln;
    sln.set_exact(mesh, source);
    this->project_global(Tuple<MeshFunction*>(&sln), Tuple<Solution*>(target), Tuple<int>(proj_norm));
  };

  /// Global orthogonal projection of one scalar ExactFunction -- user specified projection bi/linear forms.
  void project_global(ExactFunction source, Solution* target,
                  std::pair<WeakForm::biform_val_t, WeakForm::biform_ord_t> proj_biform,
                  std::pair<WeakForm::liform_val_t, WeakForm::liform_ord_t> proj_liform)
  {
    if (this->wf->neq != 1)
      error("Number of projected functions must be one if there is only one equation, in LinSystem::project_global().");
    // todo: check that supplied forms take scalar valued functions
    Mesh *mesh = this->get_mesh(0);
    if (mesh == NULL) error("Mesh is NULL in project_global().");
    Solution sln;
    sln.set_exact(mesh, source);
    this->project_global(Tuple<MeshFunction*>(&sln), Tuple<Solution*>(target), biforms_tuple_t(proj_biform), liforms_tuple_t(proj_liform));
  };

  /// Global orthogonal projection of one vector-valued ExactFunction.
  void project_global(ExactFunction2 source, Solution* target)
  {
    if (this->wf->neq != 1)
      error("Number of projected functions must be one if there is only one equation, in LinSystem::project_global().");
    int proj_norm = 2; // Hcurl
    Mesh *mesh = this->get_mesh(0);
    if (mesh == NULL) error("Mesh is NULL in project_global().");
    Solution sln;
    sln.set_exact(mesh, source);
    this->project_global(Tuple<MeshFunction*>(&sln), Tuple<Solution*>(target), Tuple<int>(proj_norm));
  };

  /// Projection-based interpolation of an exact function. This is faster than the
  /// global projection since no global matrix problem is solved.
  void project_local(ExactFunction exactfn, Mesh* mesh,
                     Solution* result, int proj_norm = H2D_DEFAULT_PROJ_NORM)
  {
    /// TODO
  }

  /// Needed for problems where BC depend on time.
  void update_essential_bc_values();

  /// Frees spaces. Called automatically on destruction.
  void free_spaces();

  Space** spaces;
  WeakForm* wf;

protected:

  CommonSolver* solver;
  CommonSolver* solver_default;

  PrecalcShapeset** pss;

  CooMatrix *A;
  bool mat_sym; ///< true if symmetric - then only upper half stored

  scalar* Vec; ///< solution coefficient vector
  int Vec_length;
  scalar* RHS; ///< assembled right-hand side
  int RHS_length;
  scalar* Dir; ///< contributions to the RHS from Dirichlet lift
  int Dir_length;

  void create_matrix(bool rhsonly);
  void insert_block(scalar** mat, int* iidx, int* jidx, int ilen, int jlen);

  ExtData<Ord>* init_ext_fns_ord(std::vector<MeshFunction *> &ext);
  ExtData<scalar>* init_ext_fns(std::vector<MeshFunction *> &ext, RefMap *rm, const int order);
  Func<double>* get_fn(PrecalcShapeset *fu, RefMap *rm, const int order);

  // Key for caching transformed function values on elements
  struct Key
  {
    int index;
    int order;
    int sub_idx;
    int shapeset_type;

    Key(int index, int order, int sub_idx, int shapeset_type)
    {
      this->index = index;
      this->order = order;
      this->sub_idx = sub_idx;
      this->shapeset_type = shapeset_type;
    }
  };

  struct Compare
  {
    bool operator()(Key a, Key b) const
    {
      if (a.index < b.index) return true;
      else if (a.index > b.index) return false;
      else
      {
        if (a.order < b.order) return true;
        else if (a.order > b.order) return false;
        else
        {
          if (a.sub_idx < b.sub_idx) return true;
          else if (a.sub_idx > b.sub_idx) return false;
          else
          {
            if (a.shapeset_type < b.shapeset_type) return true;
            else return false;
          }
        }
      }
    }
  };

  // Caching transformed values for element
  std::map<Key, Func<double>*, Compare> cache_fn;
  Geom<double>* cache_e[g_max_quad + 1 + 4];
  double* cache_jwt[g_max_quad + 1 + 4];

  void init_cache();
  void delete_cache();

  scalar eval_form(WeakForm::BiFormVol *bf, PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv);
  scalar eval_form(WeakForm::LiFormVol *lf, PrecalcShapeset *fv, RefMap *rv);
  scalar eval_form(WeakForm::BiFormSurf *bf, PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, EdgePos* ep);
  scalar eval_form(WeakForm::LiFormSurf *lf, PrecalcShapeset *fv, RefMap *rv, EdgePos* ep);

  scalar** get_matrix_buffer(int n)
  {
    if (n <= mat_size) return buffer;
    if (buffer != NULL) delete [] buffer;
    return (buffer = new_matrix<scalar>(mat_size = n));
  }

  scalar** buffer;
  int mat_size;

  int* sp_seq;
  int wf_seq;
  int num_user_pss;
  bool values_changed;
  bool struct_changed;
  bool want_dir_contrib;
  bool have_spaces;

  friend class RefSystem;

};



#endif
