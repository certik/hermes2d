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

#ifndef __HERMES2D_ADAPT_ORTHO_H1_H
#define __HERMES2D_ADAPT_ORTHO_H1_H

#include "forms.h"
#include "weakform.h"
#include "integrals_h1.h"

/// \brief hp-adaptivity module for H1 spaces.
///
/// H1OrthoHP is a fast hp-adaptivity module for continuous elements.
/// Given a reference solution, it provides functions to calculate H1 or
/// energy error estimates, acts as a container for the calculated errors
/// and contains the "ortho" hp-adaptivty algorithm based on fast
/// projections to an orthonormal set of functions.
///
class PUBLIC_API H1OrthoHP
{
public:

  /// Initializes the class. 'num' is the number of mesh-space pairs to be adapted.
  /// After 'num', exactly that many space pointers must follow.
  H1OrthoHP(int num, ...);
  virtual ~H1OrthoHP();


  typedef scalar (*biform_val_t) (int n, double *wt, Func<scalar> *u, Func<scalar> *v, Geom<double> *e, ExtData<scalar> *);
  typedef Ord (*biform_ord_t) (int n, double *wt, Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *);

  /// Sets user defined bilinear form to calculate error. Default forms are h1 error (on diagonal).
  /// Use this function only to change it (e.g. energy error).
  void set_biform(int i, int j, biform_val_t bi_form, biform_ord_t bi_ord);

  /// Type-safe version of calc_error_n() for one solution.
  virtual double calc_error(MeshFunction* sln, MeshFunction* rsln);

  /// Type-safe version of calc_error_n() for two solutions.
  double calc_error_2(MeshFunction* sln1, MeshFunction* sln2, MeshFunction* rsln1, MeshFunction* rsln2);

  /// Calculates the error of the solution using given norms. 'n' must be the
  /// same as 'num' in the constructor. After that, n coarse solution
  /// pointers are passed, followed by n fine solution pointers.
  virtual double calc_error_n(int n, ...);

  /// Selects elements to refine (based on results from calc_error() or calc_energy_error())
  /// and performs their optimal hp-refinement.
  /// \param adapt_result Contains result of adaptivity step. If NULL, result is not gathered or processed.
  bool adapt(double thr, int strat = 0, int adapt_type = 0, bool iso_only = false, int regularize = -1,
             double conv_exp = 1.0, int max_order = -1, bool same_orders = false, double to_be_processed = 0.0);

  /// Unrefines the elements with the smallest error
  void unrefine(double thr);

  /// Internal. Used by adapt(). Can be utilized in specialized adaptivity
  /// procedures, for which adapt() is not sufficient.
  /// \returns Selected candidate. If zero, no other candidate was selected.
  virtual int get_optimal_refinement(Element* e, int order, Solution* rsln, int& split, int4 p, int4 q,
                                     bool h_only = false, bool iso_only = false, double conv_exp = 1.0, 
                                     int max_order = -1);

  /// Internal. Functions to obtain errors of individual elements.
  struct ElementReference { ///< A reference to a element.
    int id, comp;
    ElementReference() {};
    ElementReference(int id, int comp) : id(id), comp(comp) {};
  };
  double get_element_error(int component, int id) const { return errors[component][id]; }
  ElementReference*  get_sorted_elements() const { return esort; }
  int    get_total_active_elements() const { return nact; }

protected: //adaptivity
  struct ElementToRefine { ///< A record of an element and a selected candidate.
    int id; //ID of element
    int comp; //componet
    int split; //proposed refinement
    int4 p; //orders
    int4 q; //H orders
    ElementToRefine() : id(-1), comp(-1) {};
    ElementToRefine(int id, int comp) : id(id), comp(comp), split(0) {};
    ElementToRefine(const ElementToRefine &orig) : id(orig.id), comp(orig.comp), split(orig.split) {
      memcpy(p, orig.p, sizeof(int4));
      memcpy(q, orig.q, sizeof(int4));
    };
    ElementToRefine& operator=(const ElementToRefine& orig) {
      id = orig.id;
      comp = orig.comp;
      split = orig.split;
      memcpy(p, orig.p, sizeof(int4));
      memcpy(q, orig.q, sizeof(int4));
      return *this;
    }
    int get_num_sons() const { ///< Returns a number of sons.
      if (split == 0)
        return 4;
      else if (split == 1 || split == 2)
        return 2;
      else
        return 1;
    };
  };

  std::queue<ElementReference> priority_esort; ///< A list of priority elements that are processed before the next element in esort is processed.

  virtual bool ignore_element_adapt(const int inx_element, const Mesh* mesh, const Element* element) { return false; }; ///< Returns true, if an element should be ignored for purposes of adaptivity.
  virtual bool can_adapt_element(Mesh* mesh, Element* e, const int split, const int4& p, const int4& q) { return true; }; ///< Returns true, if an element can be adapted using a selected candidate.
  virtual void apply_refinements(Mesh** meshes, std::vector<ElementToRefine>* elems_to_refine); ///< Apply refinements.

protected: //optimal refinement
  struct Cand ///< A candidate.
  {
    double error; ///< Error of this candidate.
    int dofs; 
    int split; ///< Operation.
    int p[4]; ///< Orders of sons.
  };

  virtual Cand* create_candidates(Element* e, int order, bool h_only, bool iso_only, int max_order, int* num_cand); ///< Creates a list of candidates. Allocated array has to be deallocated by a user of the method.
  virtual int evalute_candidates(Cand* cand, int num_cand, Element* e, int order, Solution* rsln, double* avg_error, double* dev_error); ///< Evaluates candidates. Calculates their error and dofs. Calculates average error and sample deviation.
  virtual void select_best_candidate(const Cand* cand, const int num_cand, Element* e, const double avg_error, const double dev_error, int* selected_cand, int* selected_h_cand, double conv_exp); ///< Selects the best candidate and the best h-candidate.

protected:
  // spaces & solutions
  int num;
  Space* spaces[10];
  Solution* sln[10];
  Solution* rsln[10];

  // element error arrays
  double* errors[10];
  double  norms[10]; // ?
  bool    have_errors;
  double  total_err;
  ElementReference* esort;
  int   nact;

  // bilinear forms to calculate error
  biform_val_t form[10][10];
  biform_ord_t ord[10][10];

  // evaluation of error and norm forms
  scalar eval_error(biform_val_t bi_fn, biform_ord_t bi_ord,
                    MeshFunction *sln1, MeshFunction *sln2, MeshFunction *rsln1, MeshFunction *rsln2,
                    RefMap *rv1,        RefMap *rv2,        RefMap *rrv1,        RefMap *rrv2);

  scalar eval_norm(biform_val_t bi_fn, biform_ord_t bi_ord,
                   MeshFunction *rsln1, MeshFunction *rsln2, RefMap *rrv1, RefMap *rrv2);

  // orthonormal basis tables
  static double3** obase[2][9];
  static int  basecnt[2][11];
  static bool obase_ready;

  static void calc_ortho_base();
  static int build_shape_inxs(const int mode, H1Shapeset& shapeset, int idx[121]); ///< Build indices of shape functions and initializes ranges. Returns number of indices.

  virtual void calc_projection_errors(Element* e, int order, Solution* rsln,
                                     double herr[8][11], double perr[11]); ///< Calculate various projection errors.

  /// Builds a list of elements sorted by error descending. Assumes that H1OrthoHP::errors is initialized. Initializes H1OrthoHP::esort.
  /// \param meshes: Meshes. Indices [0,.., H1OrthoHP::num-1] contains meshes of coarse solutions, indices [H1OrthoHP::num,..,2*H1OrthoHP::num - 1] contains meshes of reference solution.
  void sort_elements_by_error(Mesh** meshes);

private:
  static double** cmp_err; ///< An helper array used to sort reference to elements.
  static int compare(const void* p1, const void* p2); ///< Compares to reference to an element according to H1OrthoHP::cmp_err.

public:

  /// Internal.
  static void free_ortho_base();

};



#endif
