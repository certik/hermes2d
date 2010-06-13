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

#include "tuple.h"
#include "matrix.h"
#include "matrix_old.h"
#include "forms.h"
#include "weakform.h"
#include <map>

class Space;
class PrecalcShapeset;
class WeakForm;
class Solver;

H2D_API_USED_TEMPLATE(Tuple<Space*>); ///< Instantiated template. It is used to create a clean Windows DLL interface.
H2D_API_USED_TEMPLATE(Tuple<Solution*>); ///< Instantiated template. It is used to create a clean Windows DLL interface.
H2D_API_USED_TEMPLATE(Tuple<PrecalcShapeset*>); ///< Instantiated template. It is used to create a clean Windows DLL interface.


///
///
///
///
///
class H2D_API LinSystem
{
public:

  LinSystem();
  LinSystem(WeakForm* wf_, Solver* solver_);
  LinSystem(WeakForm* wf_);                  // solver will be set to NULL and default solver will be used
  LinSystem(WeakForm* wf_, Solver* solver_, Space* s_);
  LinSystem(WeakForm* wf_, Space* s_);       // solver will be set to NULL and default solver will be used
  LinSystem(WeakForm* wf_, Solver* solver_, Tuple<Space*> spaces_);
  LinSystem(WeakForm* wf_, Tuple<Space*> spaces_);      // solver will be set to NULL and default solver will be used
  LinSystem(WeakForm* wf_, Solver* solver_, Space* space1_, Space* space2_);

  virtual ~LinSystem();

  void init_lin(WeakForm* wf, Solver* solver);
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

  virtual void assemble(bool rhsonly = false);
  void assemble_rhs_only() { assemble(true); }
  bool solve(Tuple<Solution*> sln);
  bool solve(Solution* sln); // single equation case
  bool solve(Solution* sln1, Solution* sln2); // two equations case
  bool solve(Solution* sln1, Solution* sln2, Solution* sln3); // three equations case
  virtual void free();
  virtual void free_matrix();

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
  void get_solution_vector(std::vector<scalar>& sln_vector_out); ///< Returns a copy of the solution vector.

  /// Frees (if not NULL) and creates zero solution coefficient vectors Vec, Dir and RHS.
  void reset_coeff_vectors();

  /// Frees vectors Vec, RHS and Dir. These vectors should 
  /// not be freed anywhere else.
  void free_vectors();

  /// For debug purposes.
  void print_vector();

  /// Assigning DOF = enumerating basis functions in the FE spaces.
  int assign_dofs();  // all spaces

  /// Basic procedure performing orthogonal projection for an arbitrary number of 
  /// functions onto (the same number of) spaces determined by the LinSystem; proj_norm = 0  
  /// for L2 norm, proj_norm = 1 for H1 norm, proj_norm = 2 for Hcurl norm. Projected can 
  /// be any MeshFunction, Solution or Filter. The result of the projection will satisfy essential 
  /// boundary conditions. The projection defines the vector Vec in the class LinSystem.
  /// All projection functionality defined here is also available in the class NonlinSystem. 
  /// TODO: Implement projection-based interpolation (PBI) as an alternative of this. 
  /// PBI is almost as good as global orthogonal projection but way faster.
  void project_global(Tuple<MeshFunction*> source, Tuple<Solution*> target, int proj_norm = 1);

  /// Global orthogonal projection of MeshFunction* fn. Result of the projection is 
  /// returned as "result".
  void project_global(MeshFunction* source, Solution* target, int proj_norm = 1)
  {  project_global(Tuple<MeshFunction*>(source), Tuple<Solution*>(target), proj_norm);  }

  /// Global orthogonal projection of an exact function.
  void project_global(scalar (*exactfn)(double x, double y, scalar& dx, scalar& dy),
                      Solution* target, int proj_norm = 1)
  {
    Mesh *mesh = this->get_mesh(0);
    target->set_exact(mesh, exactfn);
    project_global(target, target, proj_norm);
  }

  /// Global orthogonal projection of two exact functions.
  void project_global(scalar (*exactfn1)(double x, double y, scalar& dx, scalar& dy), 
                      scalar (*exactfn2)(double x, double y, scalar& dx, scalar& dy),
                      Solution* target1, Solution* target2, int proj_norm = 1)
  {
    Mesh *mesh1 = this->get_mesh(0);
    Mesh *mesh2 = this->get_mesh(1);
    target1->set_exact(mesh1, exactfn1);
    target2->set_exact(mesh2, exactfn2);
    project_global(Tuple<MeshFunction*>(target1, target2), Tuple<Solution*>(target1, target2), proj_norm);
  }

  /// Global orthogonal projection of three exact functions.
  void project_global(scalar (*exactfn1)(double x, double y, scalar& dx, scalar& dy), 
                      scalar (*exactfn2)(double x, double y, scalar& dx, scalar& dy),
                      scalar (*exactfn3)(double x, double y, scalar& dx, scalar& dy),
                      Solution* target1, Solution* target2, Solution* target3, int proj_norm = 1)
  {
    Mesh *mesh1 = this->get_mesh(0);
    Mesh *mesh2 = this->get_mesh(1);
    Mesh *mesh3 = this->get_mesh(2);
    target1->set_exact(mesh1, exactfn1);
    target2->set_exact(mesh2, exactfn2);
    target3->set_exact(mesh3, exactfn3);
    project_global(Tuple<MeshFunction*>(target1, target2, target3), 
                   Tuple<Solution*>(target1, target2, target3), proj_norm);
  }

  /// Projection-based interpolation of an exact function. This is faster than the 
  /// global projection since no global matrix problem is solved. 
  void project_local(scalar (*exactfn)(double x, double y, scalar& dx, scalar& dy),
              Mesh* mesh, Solution* result, int proj_norm = 1)
  {
    /// TODO
  }

  /// Frees spaces and meshes. Called
  /// automatically on destruction.
  void free_meshes_and_spaces();

  Space** spaces;
  Mesh** meshes;
  WeakForm* wf;

protected:

  Solver* solver;
  void* slv_ctx;

  PrecalcShapeset** pss;

  CooMatrix *A;
  bool mat_sym; ///< true if symmetric and only upper half stored

  scalar* RHS; ///< assembled right-hand side
  scalar* Dir; ///< contributions to the RHS from Dirichlet DOFs
  scalar* Vec; ///< last solution vector

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
