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

#include "matrix.h"
#include "matrix_old.h"
#include "forms.h"
#include "weakform.h"
#include <map>

class Space;
class PrecalcShapeset;
class WeakForm;
class Solver;
struct Page;


///
///
///
///
///
class H2D_API LinSystem
{
public:

  LinSystem(WeakForm* wf, Solver* solver = NULL);
  virtual ~LinSystem();

  void set_spaces(int n, ...);
  void set_pss(int n, ...);
  void copy(LinSystem* sys);
  Space* get_space(int n) {
      //if (n < 0 || n >= this->wf->neq) error("Bad number of the space.");
      return this->spaces[n];
  }
  PrecalcShapeset* get_pss(int n) {
      //if (n < 0 || n >= this->wf->neq) error("Bad number of the space.");
      return this->pss[n];
  }

  void assemble(bool rhsonly = false);
  void assemble_rhs_only() { assemble(true); }
  bool solve(int n, ...);
  virtual void free();

  void save_matrix_matlab(const char* filename, const char* varname = "A");
  void save_rhs_matlab(const char* filename, const char* varname = "b");
  void save_matrix_bin(const char* filename);
  void save_rhs_bin(const char* filename);

  void enable_dir_contrib(bool enable = true) {  want_dir_contrib = enable;  }
  const scalar* get_solution_vec() const { return Vec; }

  int get_num_dofs() const { return ndofs; };
  int get_matrix_size() const;
  void get_matrix(int*& Ap, int*& Ai, scalar*& Ax, int& size) const;
  void get_rhs(scalar*& RHS, int& size) const { RHS = this->RHS; size=ndofs; }
  void get_solution_vector(scalar*& sln_vector, int& sln_vector_len) { sln_vector = Vec; sln_vector_len = ndofs; }
  void get_solution_vector(std::vector<scalar>& sln_vector_out) const; ///< Returns a copy of a solution vector.

protected:

  WeakForm* wf;
  Solver* solver;
  void* slv_ctx;

  Space** spaces;
  PrecalcShapeset** pss;

  int ndofs;
  CooMatrix *A;
  bool mat_sym; ///< true if symmetric and only upper half stored

  scalar* RHS; ///< assembled right-hand side
  scalar* Dir; ///< contributions to the RHS from Dirichlet DOFs
  scalar* Vec; ///< last solution vector

  void create_matrix(bool rhsonly);
  void precalc_sparse_structure(Page** pages);
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
