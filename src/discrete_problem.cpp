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

#include "common.h"
#include "discrete_problem.h"
#include "solver.h"
#include "traverse.h"
#include "space.h"
#include "precalc.h"
#include "shapeset_h1_all.h"
#include "refmap.h"
#include "solution.h"
#include "config.h"
#include "limit_order.h"
#include <algorithm>

#include "solvers.h"

int H2D_DEFAULT_PROJ_NORM = 1;

void qsort_int(int* pbase, size_t total_elems); // defined in qsort.cpp

//// interface /////////////////////////////////////////////////////////////////////////////////////

void DiscreteProblem::init(WeakForm* wf_, CommonSolver* solver_)
{
  if (wf_ == NULL) error("DiscreteProblem: a weak form must be given.");
  this->wf = wf_;
  this->solver_default = new CommonSolverSciPyUmfpack();
  this->solver = (solver_) ? solver_ : solver_default;
  this->wf_seq = -1;

  this->RHS = this->Dir = this->Vec = NULL;
  this->RHS_length = this->Dir_length = this->Vec_length = 0;
  this->A = NULL;
  this->mat_sym = false;

  this->spaces = NULL;
  this->pss = NULL;

  this->values_changed = true;
  this->struct_changed = true;
  this->have_spaces = false;

  this->alpha = 1.0;
}

// this is needed because of a constructor in NonlinSystem
DiscreteProblem::DiscreteProblem() {}

DiscreteProblem::DiscreteProblem(WeakForm* wf_)
{
  CommonSolver *solver_ = NULL;
  this->init(wf_, solver_);
}

DiscreteProblem::DiscreteProblem(WeakForm* wf_, Tuple<Space*> sp)
{
  CommonSolver* solver_ = NULL;
  this->init(wf_, solver_);
  this->init_spaces(sp);
  this->alloc_and_zero_vectors();
}

DiscreteProblem::DiscreteProblem(WeakForm* wf_, Space* s_)
{
  CommonSolver *solver_ = NULL;
  this->init(wf_, solver_);
  this->init_space(s_);
  this->alloc_and_zero_vectors();
}

DiscreteProblem::~DiscreteProblem()
{
  /* FIXME - this should be uncommented but then it gives double-free
             segfaults in adaptive examples.
  free();
  if (this->sp_seq != NULL) delete [] this->sp_seq;
  if (this->pss != NULL) delete [] this->pss;
  if (this->solver != NULL) this->solver->free_context(this->slv_ctx);
  */
  free_vectors();
  delete this->solver_default;
}

void DiscreteProblem::free_spaces()
{
  // free spaces, making sure that duplicated ones do not get deleted twice
  if (this->spaces != NULL)
  {
    for (int i = 0; i < this->wf->neq; i++) {
      // this loop skipped if there is only one space
      for (int j = i+1; j < this->wf->neq; j++) {
        if (this->spaces[j] == this->spaces[i]) this->spaces[j] = NULL;
      }
    }
    for (int i = 0; i < this->wf->neq; i++) {
      if (this->spaces[i] != NULL) {
        delete this->spaces[i];
        this->spaces[i] = NULL;
      }
    }
    delete this->spaces;
    this->spaces = NULL;
  }
}

// Should not be called by the user.
void DiscreteProblem::init_spaces(Tuple<Space*> sp)
{
  int n = sp.size();
  if (n != this->wf->neq)
    error("Number of spaces does not match number of equations in DiscreteProblem::init_spaces().");

  // initialize spaces
  this->spaces = new Space*[this->wf->neq];
  for (int i = 0; i < this->wf->neq; i++) this->spaces[i] = sp[i];
  this->sp_seq = new int[this->wf->neq];
  memset(sp_seq, -1, sizeof(int) * this->wf->neq);
  this->assign_dofs(); // Create global enumeration of DOF in all spacesin the system. NOTE: this
                       // overwrites possible existing local enumeration of DOF in the spaces
  this->have_spaces = true;

  // initialize precalc shapesets
  this->pss = new PrecalcShapeset*[this->wf->neq];
  for (int i=0; i<this->wf->neq; i++) this->pss[i] = NULL;
  this->num_user_pss = 0;
  for (int i = 0; i < n; i++){
    Shapeset *shapeset = spaces[i]->get_shapeset();
    if (shapeset == NULL) error("Internal in DiscreteProblem::init_spaces().");
    PrecalcShapeset *p = new PrecalcShapeset(shapeset);
    if (p == NULL) error("New PrecalcShapeset could not be allocated in DiscreteProblem::init_spaces().");
    this-> pss[i] = p;
    this->num_user_pss++;
  }
}

// Should not be called by the user.
void DiscreteProblem::init_space(Space* s)
{
  if (this->wf->neq != 1)
    error("Do not call init_space() for PDE systems, call init_spaces() instead.");
  this->init_spaces(Tuple<Space*>(s));
}

// Obsolete. Should be removed after FeProblem is removed.
void DiscreteProblem::set_spaces(Tuple<Space*>spaces)
{
  this->init_spaces(spaces);
}

void DiscreteProblem::set_pss(Tuple<PrecalcShapeset*> pss)
{
  warn("Call to deprecated function DiscreteProblem::set_pss().");
  int n = pss.size();
  if (n != this->wf->neq)
    error("The number of precalculated shapesets must match the number of equations.");

  for (int i = 0; i < n; i++) this->pss[i] = pss[i];
  num_user_pss = n;
}

void DiscreteProblem::set_pss(PrecalcShapeset* pss)
{
  this->set_pss(Tuple<PrecalcShapeset*>(pss));
}

void DiscreteProblem::copy(DiscreteProblem* sys)
{
  error("Not implemented yet.");
}

void DiscreteProblem::free_vectors()
{
  if (this->Vec != NULL || this->RHS != NULL || this->Dir != NULL)
    //printf("debug: freeing vectors Vec, RHS, Dir for lengths   %d\n", this->Vec_length);

  if (this->RHS != NULL) {
    delete [] this->RHS;
    this->RHS = NULL;
    this->RHS_length = 0;
  }
  if (this->Dir != NULL) {
    delete [] this->Dir;
    this->Dir = NULL;
    this->Dir_length = 0;
  }
  if (this->Vec != NULL) {
    delete [] this->Vec;
    this->Vec = NULL;
    this->Vec_length = 0;
  }
}

void DiscreteProblem::alloc_and_zero_vectors()
{
  int ndof = this->get_num_dofs();
  //printf("debug: allocating vectors Vec, RHS, Dir for ndof   %d\n", ndof);

  if (this->RHS != NULL || this->Dir != NULL || this->Vec != NULL)
    error("All vectors must be NULL in alloc_and_zero_vectors() to prevent loss of information.");

  this->Vec = new scalar[ndof];
  if (Vec == NULL) error("Not enough memory DiscreteProblem::alloc_and_zero_vectors().");
  memset(this->Vec, 0, ndof*sizeof(scalar));
  this->Vec_length = ndof;

  this->RHS = new scalar[ndof];
  if (RHS == NULL) error("Not enough memory in DiscreteProblem::alloc_and_zero_vectors().");
  memset(this->RHS, 0, ndof*sizeof(scalar));
  this->RHS_length = ndof;

  this->Dir = new scalar[ndof];
  if (Dir == NULL) error("Not enough memory in DiscreteProblem::alloc_and_zero_vectors().");
  memset(this->Dir, 0, ndof*sizeof(scalar));
  this->Dir_length = ndof;
}

void DiscreteProblem::realloc_and_zero_vectors()
{
  int ndof = this->get_num_dofs();
  //printf("debug: reallocating vectors Vec, RHS, Dir length   %d -> %d\n", this->Vec_length, ndof);

  this->Vec = (scalar*)realloc(this->Vec, ndof*sizeof(scalar));
  if (this->Vec == NULL) error("Not enough memory DiscreteProblem::realloc_and_zero_vectors().");
  memset(this->Vec, 0, ndof*sizeof(scalar));
  this->Vec_length = ndof;

  this->RHS = (scalar*)realloc(this->RHS, ndof*sizeof(scalar));
  if (this->RHS == NULL) error("Not enough memory DiscreteProblem::realloc_and_zero_vectors().");
  memset(this->RHS, 0, ndof*sizeof(scalar));
  this->RHS_length = ndof;

  this->Dir = (scalar*)realloc(this->Dir, ndof*sizeof(scalar));
  if (this->Dir == NULL) error("Not enough memory DiscreteProblem::realloc_and_zero_vectors().");
  memset(this->Dir, 0, ndof*sizeof(scalar));
  this->Dir_length = ndof;
}

void DiscreteProblem::free()
{
  free_matrix();
  free_vectors();
  free_spaces();

  this->struct_changed = this->values_changed = true;
  memset(this->sp_seq, -1, sizeof(int) * this->wf->neq);
  this->wf_seq = -1;
}

void DiscreteProblem::free_matrix()
{
  if (this->A != NULL) { ::delete this->A; this->A = NULL; }
}

//// matrix creation ///////////////////////////////////////////////////////////////////////////////

void DiscreteProblem::create_matrix(bool rhsonly)
{
  // sanity checks
  if (this->wf == NULL) error("this->wf is NULL in DiscreteProblem::get_num_dofs().");
  if (this->wf->neq == 0) error("this->wf->neq is 0 in DiscreteProblem::get_num_dofs().");
  if (this->spaces == NULL) error("this->spaces[%d] is NULL in DiscreteProblem::get_num_dofs().");

  // check if we can reuse the matrix structure
  bool up_to_date = true;
  for (int i = 0; i < this->wf->neq; i++) {
    if (this->spaces[i] ==  NULL) error("this->spaces[%d] is NULL in DiscreteProblem::get_num_dofs().", i);
    if (this->spaces[i]->get_seq() != this->sp_seq[i]) {
      up_to_date = false;
      break;
    }
  }
  if (this->wf->get_seq() != this->wf_seq) up_to_date = false;

  // calculate the number of DOF
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof = 0 in DiscreteProblem::create_matrix().");

  // if the matrix has not changed, just zero the values and we're done
  if (up_to_date)
  {
    verbose("Reusing matrix sparse structure.");
    if (!rhsonly) {
#ifdef H2D_COMPLEX
      this->A = new CooMatrix(ndof, true);
#else
      this->A = new CooMatrix(ndof);
#endif
      memset(this->Dir, 0, sizeof(scalar) * ndof);
    }
    memset(this->RHS, 0, sizeof(scalar) * ndof);
    return;
  }
  else if (rhsonly)
    error("Cannot reassemble RHS only: spaces have changed.");

  // spaces have changed: create the matrix from scratch
  this->free_matrix();
  trace("Creating matrix sparse structure...");
  TimePeriod cpu_time;

#ifdef H2D_COMPLEX
  this->A = new CooMatrix(ndof, true);
#else
  this->A = new CooMatrix(ndof);
#endif

  // save space seq numbers and weakform seq number, so we can detect their changes
  for (int i = 0; i < this->wf->neq; i++)
    this->sp_seq[i] = this->spaces[i]->get_seq();
  this->wf_seq = this->wf->get_seq();

  this->struct_changed = true;
}


int DiscreteProblem::get_matrix_size()
{
    return this->A->get_size();
}

void DiscreteProblem::get_matrix(int*& Ap, int*& Ai, scalar*& Ax, int& size)
{
    /// FIXME: this is a memory leak:
    CSRMatrix *m = new CSRMatrix(this->A);
    Ap = m->get_Ap(); Ai = m->get_Ai();
    size = m->get_size();
#ifdef H2D_COMPLEX
    Ax = m->get_Ax_cplx();
#else
    Ax = m->get_Ax();
#endif
}


//// assembly //////////////////////////////////////////////////////////////////////////////////////

void DiscreteProblem::insert_block(Matrix *A, scalar** mat, int* iidx, int* jidx, int ilen, int jlen)
{
    A->add_block(iidx, ilen, jidx, jlen, mat);
}

void DiscreteProblem::assemble(Matrix* &mat_ext, scalar* &dir_ext, scalar* &rhs_ext, bool rhsonly)
{
  // sanity checks
  if (this->have_spaces == false)
    error("Before assemble(), you need to initialize spaces.");
  if (this->wf == NULL) error("this->wf = NULL in DiscreteProblem::assemble().");
  if (this->spaces == NULL) error("this->spaces = NULL in DiscreteProblem::assemble().");
  int n = this->wf->neq;
  for (int i=0; i<n; i++) if (this->spaces[i] == NULL)
			    error("this->spaces[%d] is NULL in DiscreteProblem::assemble().", i);

  // enumerate DOF to get new length of the vectors Vec, RHS and Dir,
  // and realloc these vectors if needed
  this->assign_dofs();
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof = 0 in DiscreteProblem::assemble().");
  if (this->Vec_length != ndof || this->RHS_length != ndof || this->Dir_length != ndof) {
    this->realloc_and_zero_vectors();
  }

  if (!rhsonly) free_matrix();
  int k, m, marker;
  std::vector<AsmList> al(wf->neq);
  AsmList* am, * an;
  bool bnd[4];
  std::vector<bool> nat(wf->neq), isempty(wf->neq);
  EdgePos ep[4];
  reset_warn_order();

  if (rhsonly && this->A == NULL)
    error("Cannot reassemble RHS only: matrix is has not been assembled yet.");

  // create the sparse structure
  create_matrix(rhsonly);

  trace("Assembling stiffness matrix...");
  TimePeriod cpu_time;

  // create slave pss's for test functions, init quadrature points
  std::vector<PrecalcShapeset*> spss(wf->neq, static_cast<PrecalcShapeset*>(NULL));
  PrecalcShapeset *fu, *fv;
  std::vector<RefMap> refmap(wf->neq);
  for (int i = 0; i < wf->neq; i++)
  {
    spss[i] = new PrecalcShapeset(pss[i]);
    pss [i]->set_quad_2d(&g_quad_2d_std);
    spss[i]->set_quad_2d(&g_quad_2d_std);
    refmap[i].set_quad_2d(&g_quad_2d_std);
  }

  // initialize buffer
  buffer = NULL;
  mat_size = 0;
  get_matrix_buffer(9);

  // obtain a list of assembling stages
  std::vector<WeakForm::Stage> stages;
  wf->get_stages(spaces, stages, rhsonly);

  // Loop through all assembling stages -- the purpose of this is increased performance
  // in multi-mesh calculations, where, e.g., only the right hand side uses two meshes.
  // In such a case, the matrix forms are assembled over one mesh, and only the rhs
  // traverses through the union mesh. On the other hand, if you don't use multi-mesh
  // at all, there will always be only one stage in which all forms are assembled as usual.
  Traverse trav;
  for (unsigned int ss = 0; ss < stages.size(); ss++)
  {
    WeakForm::Stage* s = &stages[ss];
    for (unsigned int i = 0; i < s->idx.size(); i++)
      s->fns[i] = pss[s->idx[i]];
    for (unsigned int i = 0; i < s->ext.size(); i++)
      s->ext[i]->set_quad_2d(&g_quad_2d_std);
    trav.begin(s->meshes.size(), &(s->meshes.front()), &(s->fns.front()));

    // assemble one stage
    Element** e;
    while ((e = trav.get_next_state(bnd, ep)) != NULL)
    {
      // find a non-NULL e[i]
      Element* e0 = NULL;
      for (unsigned int i = 0; i < s->idx.size(); i++)
        if ((e0 = e[i]) != NULL) break;
      if (e0 == NULL) continue;

      // set maximum integration order for use in integrals, see limit_order()
      update_limit_table(e0->get_mode());

      // obtain assembly lists for the element at all spaces, set appropriate mode for each pss
      std::fill(isempty.begin(), isempty.end(), false);
      for (unsigned int i = 0; i < s->idx.size(); i++)
      {
        int j = s->idx[i];
        if (e[i] == NULL) { isempty[j] = true; continue; }
        spaces[j]->get_element_assembly_list(e[i], &al[j]);
        /** \todo Do not retrieve assembly list again if the element has not changed */

        spss[j]->set_active_element(e[i]);
        spss[j]->set_master_transform();
        refmap[j].set_active_element(e[i]);
        refmap[j].force_transform(pss[j]->get_transform(), pss[j]->get_ctm());
      }
      marker = e0->marker;

      init_cache();
      //// assemble volume matrix forms //////////////////////////////////////
      for (unsigned int ww = 0; ww < s->mfvol.size(); ww++)
      {
        WeakForm::MatrixFormVol* mfv = s->mfvol[ww];
        if (isempty[mfv->i] || isempty[mfv->j]) continue;
        if (mfv->area != H2D_ANY && !wf->is_in_area(marker, mfv->area)) continue;
        m = mfv->i;  fv = spss[m];  am = &al[m];
        n = mfv->j;  fu = pss[n];   an = &al[n];
        bool tra = (m != n) && (mfv->sym != 0);
        bool sym = (m == n) && (mfv->sym == 1);

        // assemble the local stiffness matrix for the form mfv
        scalar bi, **mat = get_matrix_buffer(std::max(am->cnt, an->cnt));
        for (int i = 0; i < am->cnt; i++)
        {
          if (!tra && (k = am->dof[i]) < 0) continue;
          fv->set_active_shape(am->idx[i]);

          if (!sym) // unsymmetric block
          {
            for (int j = 0; j < an->cnt; j++) {
              fu->set_active_shape(an->idx[j]);
              // FIXME - the NULL on the following line is temporary, an array of solutions 
              // should be passed there.
              bi = eval_form(mfv, NULL, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
              if (an->dof[j] < 0) Dir[k] -= bi; else mat[i][j] = bi;
            }
          }
          else // symmetric block
          {
            for (int j = 0; j < an->cnt; j++) {
              if (j < i && an->dof[j] >= 0) continue;
              fu->set_active_shape(an->idx[j]);
              // FIXME - the NULL on the following line is temporary, an array of solutions 
              // should be passed there.
              bi = eval_form(mfv, NULL, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
              if (an->dof[j] < 0) Dir[k] -= bi; else mat[i][j] = mat[j][i] = bi;
            }
          }
        }

        // insert the local stiffness matrix into the global one
        insert_block(this->A, mat, am->dof, an->dof, am->cnt, an->cnt);

        // insert also the off-diagonal (anti-)symmetric block, if required
        if (tra)
        {
          if (mfv->sym < 0) chsgn(mat, am->cnt, an->cnt);
          transpose(mat, am->cnt, an->cnt);
          insert_block(this->A, mat, an->dof, am->dof, an->cnt, am->cnt);

          // we also need to take care of the RHS...
          for (int j = 0; j < am->cnt; j++)
            if (am->dof[j] < 0)
              for (int i = 0; i < an->cnt; i++)
                if (an->dof[i] >= 0)
                  Dir[an->dof[i]] -= mat[i][j];
        }
      }

      //// assemble volume linear forms ////////////////////////////////////////
      for (unsigned int ww = 0; ww < s->vfvol.size(); ww++)
      {
        WeakForm::VectorFormVol* vfv = s->vfvol[ww];
        if (isempty[vfv->i]) continue;
        if (vfv->area != H2D_ANY && !wf->is_in_area(marker, vfv->area)) continue;
        m = vfv->i;  fv = spss[m];  am = &al[m];

        for (int i = 0; i < am->cnt; i++)
        {
          if (am->dof[i] < 0) continue;
          fv->set_active_shape(am->idx[i]);
          // FIXME - the NULL on the following line is temporary, an array of solutions 
          // should be passed there.
          RHS[am->dof[i]] += eval_form(vfv, NULL, fv, &refmap[m]) * am->coef[i];
        }
      }


      // assemble surface integrals now: loop through boundary edges of the element
      for (unsigned int edge = 0; edge < e0->nvert; edge++)
      {
        if (!bnd[edge]) continue;
        marker = ep[edge].marker;

        // obtain the list of shape functions which are nonzero on this edge
        for (unsigned int i = 0; i < s->idx.size(); i++) {
          if (e[i] == NULL) continue;
          int j = s->idx[i];
          if ((nat[j] = (spaces[j]->bc_type_callback(marker) == BC_NATURAL)))
            spaces[j]->get_edge_assembly_list(e[i], edge, &al[j]);
        }

        // assemble surface matrix forms ///////////////////////////////////
        for (unsigned int ww = 0; ww < s->mfsurf.size(); ww++)
        {
          WeakForm::MatrixFormSurf* mfs = s->mfsurf[ww];
          if (isempty[mfs->i] || isempty[mfs->j]) continue;
          if (mfs->area != H2D_ANY && !wf->is_in_area(marker, mfs->area)) continue;
          m = mfs->i;  fv = spss[m];  am = &al[m];
          n = mfs->j;  fu = pss[n];   an = &al[n];

          if (!nat[m] || !nat[n]) continue;
          ep[edge].base = trav.get_base();
          ep[edge].space_v = spaces[m];
          ep[edge].space_u = spaces[n];

          scalar bi, **mat = get_matrix_buffer(std::max(am->cnt, an->cnt));
          for (int i = 0; i < am->cnt; i++)
          {
            if ((k = am->dof[i]) < 0) continue;
            fv->set_active_shape(am->idx[i]);
            for (int j = 0; j < an->cnt; j++)
            {
              fu->set_active_shape(an->idx[j]);
              // FIXME - the NULL on the following line is temporary, an array of solutions 
              // should be passed there.
              bi = eval_form(mfs, NULL, fu, fv, &refmap[n], &refmap[m], &(ep[edge])) * an->coef[j] * am->coef[i];
              if (an->dof[j] >= 0) mat[i][j] = bi; else Dir[k] -= bi;
            }
          }
          insert_block(this->A, mat, am->dof, an->dof, am->cnt, an->cnt);
        }

        // assemble surface linear forms /////////////////////////////////////
        for (unsigned int ww = 0; ww < s->vfsurf.size(); ww++)
        {
          WeakForm::VectorFormSurf* vfs = s->vfsurf[ww];
          if (isempty[vfs->i]) continue;
          if (vfs->area != H2D_ANY && !wf->is_in_area(marker, vfs->area)) continue;
          m = vfs->i;  fv = spss[m];  am = &al[m];

          if (!nat[m]) continue;
          ep[edge].base = trav.get_base();
          ep[edge].space_v = spaces[m];

          for (int i = 0; i < am->cnt; i++)
          {
            if (am->dof[i] < 0) continue;
            fv->set_active_shape(am->idx[i]);
            // FIXME - the NULL on the following line is temporary, an array of solutions 
            // should be passed there.
            RHS[am->dof[i]] += eval_form(vfs, NULL, fv, &refmap[m], &(ep[edge])) * am->coef[i];
          }
        }
      }
      delete_cache();
    }
    trav.finish();
  }

  verbose("Stiffness matrix assembled (stages: %d)", stages.size());
  report_time("Stiffness matrix assembled in %g s", cpu_time.tick().last());
  for (int i = 0; i < wf->neq; i++) delete spss[i];
  delete [] buffer;

  if (!rhsonly) values_changed = true;

  mat_ext = this->A;
  dir_ext = this->Dir;
  rhs_ext = this->RHS;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Initialize integration order for external functions
ExtData<Ord>* DiscreteProblem::init_ext_fns_ord(std::vector<MeshFunction *> &ext)
{
  ExtData<Ord>* fake_ext = new ExtData<Ord>;
  fake_ext->nf = ext.size();
  Func<Ord>** fake_ext_fn = new Func<Ord>*[fake_ext->nf];
  for (int i = 0; i < fake_ext->nf; i++)
    fake_ext_fn[i] = init_fn_ord(ext[i]->get_fn_order());
  fake_ext->fn = fake_ext_fn;

  return fake_ext;
}

// Initialize external functions (obtain values, derivatives,...)
ExtData<scalar>* DiscreteProblem::init_ext_fns(std::vector<MeshFunction *> &ext, RefMap *rm, const int order)
{
  ExtData<scalar>* ext_data = new ExtData<scalar>;
  Func<scalar>** ext_fn = new Func<scalar>*[ext.size()];
  for (unsigned int i = 0; i < ext.size(); i++)
    ext_fn[i] = init_fn(ext[i], rm, order);
  ext_data->nf = ext.size();
  ext_data->fn = ext_fn;

  return ext_data;

}

// Initialize shape function values and derivatives (fill in the cache)
Func<double>* DiscreteProblem::get_fn(PrecalcShapeset *fu, RefMap *rm, const int order)
{
  Key key(256 - fu->get_active_shape(), order, fu->get_transform(), fu->get_shapeset()->get_id());
  if (cache_fn[key] == NULL)
    cache_fn[key] = init_fn(fu, rm, order);

  return cache_fn[key];
}

// Caching transformed values
void DiscreteProblem::init_cache()
{
  for (int i = 0; i < g_max_quad + 1 + 4; i++)
  {
    cache_e[i] = NULL;
    cache_jwt[i] = NULL;
  }
}

void DiscreteProblem::delete_cache()
{
  for (int i = 0; i < g_max_quad + 1 + 4; i++)
  {
    if (cache_e[i] != NULL)
    {
      cache_e[i]->free(); delete cache_e[i];
      delete [] cache_jwt[i];
    }
  }
  for (std::map<Key, Func<double>*, Compare>::iterator it = cache_fn.begin(); it != cache_fn.end(); it++)
  {
    (it->second)->free_fn(); delete (it->second);
  }
  cache_fn.clear();
}

//// evaluation of forms, general case ///////////////////////////////////////////////////////////


// Actual evaluation of volume Jacobian form (calculates integral)
scalar DiscreteProblem::eval_form(WeakForm::MatrixFormVol *mfv, Solution *sln[], 
                        PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv)
{
  // determine the integration order
  int inc = (fu->get_num_components() == 2) ? 1 : 0;
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  //for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
  if (sln != NULL) {
    for (int i = 0; i < wf->neq; i++) {
      if (sln[i] != NULL) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
      else oi[i] = init_fn_ord(0);
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(0);
  }

  Func<Ord>* ou = init_fn_ord(fu->get_fn_order() + inc);
  Func<Ord>* ov = init_fn_ord(fv->get_fn_order() + inc);
  ExtData<Ord>* fake_ext = init_ext_fns_ord(mfv->ext);

  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  Ord o = mfv->ord(1, &fake_wt, oi, ou, ov, fake_e, fake_ext);
  int order = ru->get_inv_ref_order();
  order += o.get_order();
  limit_order_nowarn(order);

  for (int i = 0; i < wf->neq; i++) {  
    if (oi[i] != NULL) { oi[i]->free_ord(); delete oi[i]; }
  }
  if (ou != NULL) {
    ou->free_ord(); delete ou;
  }
  if (ov != NULL) {
    ov->free_ord(); delete ov;
  }
  if (fake_e != NULL) delete fake_e;
  if (fake_ext != NULL) {fake_ext->free_ord(); delete fake_ext;}

  // eval the form
  Quad2D* quad = fu->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // init geometry and jacobian*weights
  if (cache_e[order] == NULL)
  {
    cache_e[order] = init_geom_vol(ru, order);
    double* jac = ru->get_jacobian(order);
    cache_jwt[order] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[order][i] = pt[i][2] * jac[i];
  }
  Geom<double>* e = cache_e[order];
  double* jwt = cache_jwt[order];

  // function values and values of external functions
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  //for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(sln[i], rv, order);
  if (sln != NULL) {
    for (int i = 0; i < wf->neq; i++) {
      if (sln[i] != NULL) prev[i]  = init_fn(sln[i], rv, order);
      else prev[i] = NULL;
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) prev[i] = NULL;
  }

  Func<double>* u = get_fn(fu, ru, order);
  Func<double>* v = get_fn(fv, rv, order);
  ExtData<scalar>* ext = init_ext_fns(mfv->ext, rv, order);

  scalar res = mfv->fn(np, jwt, prev, u, v, e, ext);

  for (int i = 0; i < wf->neq; i++) {  
    if (prev[i] != NULL) prev[i]->free_fn(); delete prev[i]; 
  }
  if (ext != NULL) {ext->free(); delete ext;}
  return res;
}


// Actual evaluation of volume vector form (calculates integral)
scalar DiscreteProblem::eval_form(WeakForm::VectorFormVol *vfv, Solution *sln[], PrecalcShapeset *fv, RefMap *rv)
{
  // determine the integration order
  int inc = (fv->get_num_components() == 2) ? 1 : 0;
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  //for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
  if (sln != NULL) {
    for (int i = 0; i < wf->neq; i++) {
      if (sln[i] != NULL) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
      else oi[i] = init_fn_ord(0);
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(0);
  }
  Func<Ord>* ov = init_fn_ord(fv->get_fn_order() + inc);
  ExtData<Ord>* fake_ext = init_ext_fns_ord(vfv->ext);

  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  Ord o = vfv->ord(1, &fake_wt, oi, ov, fake_e, fake_ext);
  int order = rv->get_inv_ref_order();
  order += o.get_order();
  limit_order_nowarn(order);

  for (int i = 0; i < wf->neq; i++) { 
    if (oi[i] != NULL) {
      oi[i]->free_ord(); delete oi[i]; 
    }
  }
  if (ov != NULL) {ov->free_ord(); delete ov;}
  if (fake_e != NULL) delete fake_e;
  if (fake_ext != NULL) {fake_ext->free_ord(); delete fake_ext;}

  // eval the form
  Quad2D* quad = fv->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // init geometry and jacobian*weights
  if (cache_e[order] == NULL)
  {
    cache_e[order] = init_geom_vol(rv, order);
    double* jac = rv->get_jacobian(order);
    cache_jwt[order] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[order][i] = pt[i][2] * jac[i];
  }
  Geom<double>* e = cache_e[order];
  double* jwt = cache_jwt[order];

  // function values and values of external functions
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  //for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(sln[i], rv, order);
  if (sln != NULL) {
    for (int i = 0; i < wf->neq; i++) {
      if (sln[i] != NULL) prev[i]  = init_fn(sln[i], rv, order);
      else prev[i] = NULL;
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) prev[i] = NULL;
  }

  Func<double>* v = get_fn(fv, rv, order);
  ExtData<scalar>* ext = init_ext_fns(vfv->ext, rv, order);

  scalar res = vfv->fn(np, jwt, prev, v, e, ext);

  for (int i = 0; i < wf->neq; i++) { 
    if (prev[i] != NULL) {
      prev[i]->free_fn(); delete prev[i]; 
    }
  }
  if (ext != NULL) {ext->free(); delete ext;}
  return res;

}

// Actual evaluation of surface Jacobian form (calculates integral)
scalar DiscreteProblem::eval_form(WeakForm::MatrixFormSurf *mfs, Solution *sln[], 
                        PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, EdgePos* ep)
{
  // eval the form
  Quad2D* quad = fu->get_quad_2d();
  // FIXME - this needs to be order-dependent
  int eo = quad->get_edge_points(ep->edge);
  double3* pt = quad->get_points(eo);
  int np = quad->get_num_points(eo);

  // init geometry and jacobian*weights
  if (cache_e[eo] == NULL)
  {
    cache_e[eo] = init_geom_surf(ru, ep, eo);
    double3* tan = ru->get_tangent(ep->edge);
    cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }
  Geom<double>* e = cache_e[eo];
  double* jwt = cache_jwt[eo];

  // function values and values of external functions
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  //for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(sln[i], rv, eo);
  if (sln != NULL) {
    for (int i = 0; i < wf->neq; i++) {
      if (sln[i] != NULL) prev[i]  = init_fn(sln[i], rv, eo);
      else prev[i] = NULL;
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) prev[i] = NULL;
  }

  Func<double>* u = get_fn(fu, ru, eo);
  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(mfs->ext, rv, eo);

  scalar res = mfs->fn(np, jwt, prev, u, v, e, ext);

  for (int i = 0; i < wf->neq; i++) { 
    if (prev[i] != NULL) {
      prev[i]->free_fn(); delete prev[i]; 
    }
  }
  if (ext != NULL) {ext->free(); delete ext;}
  return 0.5 * res; // Edges are parameterized from 0 to 1 while integration weights
                    // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                    // the weights.
}


// Actual evaluation of surface vector form (calculates integral)
scalar DiscreteProblem::eval_form(WeakForm::VectorFormSurf *vfs, Solution *sln[], 
                        PrecalcShapeset *fv, RefMap *rv, EdgePos* ep)
{
  // eval the form
  Quad2D* quad = fv->get_quad_2d();
  // FIXME - this needs to be order-dependent
  int eo = quad->get_edge_points(ep->edge);
  double3* pt = quad->get_points(eo);
  int np = quad->get_num_points(eo);

  // init geometry and jacobian*weights
  if (cache_e[eo] == NULL)
  {
    cache_e[eo] = init_geom_surf(rv, ep, eo);
    double3* tan = rv->get_tangent(ep->edge);
    cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }
  Geom<double>* e = cache_e[eo];
  double* jwt = cache_jwt[eo];

  // function values and values of external functions
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  //for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(sln[i], rv, eo);
  if (sln != NULL) {
    for (int i = 0; i < wf->neq; i++) {
      if (sln[i] != NULL) prev[i]  = init_fn(sln[i], rv, eo);
      else prev[i] = NULL;
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) prev[i] = NULL;
  }

  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(vfs->ext, rv, eo);

  scalar res = vfs->fn(np, jwt, prev, v, e, ext);

  for (int i = 0; i < wf->neq; i++) {  
    if (prev[i] != NULL) {prev[i]->free_fn(); delete prev[i]; }
  }
  if (ext != NULL) {ext->free(); delete ext;}
  return 0.5 * res; // Edges are parameterized from 0 to 1 while integration weights
                    // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                    // the weights.
}



//// solve /////////////////////////////////////////////////////////////////////////////////////////

bool DiscreteProblem::solve_matrix_problem(Matrix* mat, scalar* vec) 
{
  // check matrix size
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof = 0 in DiscreteProblem::solve().");
  if (ndof != this->A->get_size())
    error("Matrix size does not match vector length in DiscreteProblem:solve().");

  // FIXME: similar test should be done for the vector "vec" also, but we need
  // to access the information about its length.
  // ...

  // solve the matrix problem (and report time)
  TimePeriod cpu_time;
  bool flag = this->solver->solve(mat, vec);
  report_time("Matrix problem solved in %g s", cpu_time.tick().last());

  return flag;
}

bool DiscreteProblem::solve(Matrix* mat, scalar* rhs, scalar* vec)
{
  int ndof = this->get_num_dofs();

  // sanity checks
  if (mat == NULL) error("matrix is NULL in DiscreteProblem::solve().");
  if (rhs == NULL) error("rhs is NULL in DiscreteProblem::solve().");
  if (vec == NULL) error("vec is NULL in DiscreteProblem::solve().");
  if (ndof == 0) error("ndof = 0 in DiscreteProblem::solve().");
  if (ndof != mat->get_size())
    error("Matrix size does not match ndof in in DiscreteProblem:solve().");

  // copy "vec" into "delta" and solve the matrix problem with "mat", "delta"
  scalar* delta = new scalar[ndof];
  memcpy(delta, rhs, sizeof(scalar) * ndof);
  bool flag = this->solve_matrix_problem(mat, delta);
  if (flag == false) return false;

  // add the result which is in "delta" to the previous 
  // solution vector which is in "vec"
  for (int i = 0; i < ndof; i++) vec[i] += delta[i];
  delete [] delta;

  return true;
}

bool DiscreteProblem::solve(Tuple<Solution*> sln)
{
  int ndof = this->get_num_dofs();
  int n = sln.size();
 
  // sanity checks
  if (n != this->wf->neq)
    error("Number of solutions does not match the number of equations in DiscreteProblem::solve().");
  if (this->Vec == NULL) error("Vec is NULL in DiscreteProblem::solve().");
  if (this->Vec_length != ndof || this->RHS_length != ndof || this->Dir_length != ndof)
    error("Length of vectors Vec, RHS or Dir does not match this->ndof in DiscreteProblem::solve().");
  if (ndof == 0) error("ndof = 0 in DiscreteProblem::solve().");
  if (ndof != this->A->get_size())
    error("Matrix size does not match vector length in DiscreteProblem:solve().");

  // solve the matrix problem with this->A and this->RHS, and add the
  // result to this->Vec
  bool flag = this->solve(this->A, this->RHS, this->Vec);
  if (flag == false) return false; 

  // copy this->Vec into Solutions
  if (this->spaces == NULL) error("this->spaces == NULL in DiscreteProblem::solve().");
  if (this->pss == NULL) error("this->pss == NULL in DiscreteProblem::solve().");
  if (this->Vec == NULL) error("this->Vec == NULL in LinearProblem::solve().");
  for (int i = 0; i < n; i++)
  {
    if(this->spaces[i] == NULL) error("this->spaces[%d] == NULL in LinearProblem::solve().", i);
    if(this->spaces[i]->get_mesh() == NULL) error("this->spaces[%d]->get_mesh() == NULL in LinearProblem::solve().", i);
    sln[i]->set_fe_solution(this->spaces[i], this->pss[i], this->Vec);
  }

  return true;
}

//// matrix and solution output /////////////////////////////////////////////////////////////////////////////////

void DiscreteProblem::get_solution_vector(std::vector<scalar>& sln_vector_out) {
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof = 0 in DiscreteProblem::get_solution_vector().");

  std::vector<scalar> temp(ndof);
  sln_vector_out.swap(temp);
  for(int i = 0; i < ndof; i++)
    sln_vector_out[i] = Vec[i];
}

void DiscreteProblem::save_matrix_matlab(const char* filename, const char* varname)
{
  warn("Saving matrix in Matlab format not implemented yet.");
}

void DiscreteProblem::save_rhs_matlab(const char* filename, const char* varname)
{
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof = 0 in DiscreteProblem::save_rhs_matlab().");
  if (RHS == NULL) error("RHS has not been assembled yet.");
  FILE* f = fopen(filename, "w");
  if (f == NULL) error("Could not open file %s for writing.", filename);
  verbose("Saving RHS vector in MATLAB format...");
  fprintf(f, "%% Size: %dx1\n%s = [\n", ndof, varname);
  for (int i = 0; i < ndof; i++)
    #ifndef H2D_COMPLEX
      fprintf(f, "%.18e\n", RHS[i]);
    #else
      fprintf(f, "%.18e + %.18ei\n", RHS[i].real(), RHS[i].imag());
    #endif
  fprintf(f, "];\n");
  fclose(f);
}


void DiscreteProblem::save_matrix_bin(const char* filename)
{
}

void DiscreteProblem::save_rhs_bin(const char* filename)
{
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof = 0 in DiscreteProblem::save_rhs_bin().");
  if (RHS == NULL) error("RHS has not been assembled yet.");
  FILE* f = fopen(filename, "wb");
  if (f == NULL) error("Could not open file %s for writing.", filename);
  verbose("Saving RHS vector in binary format...");
  hermes2d_fwrite("H2DR\001\000\000\000", 1, 8, f);
  int ssize = sizeof(scalar);
  hermes2d_fwrite(&ssize, sizeof(int), 1, f);
  hermes2d_fwrite(&ndof, sizeof(int), 1, f);
  hermes2d_fwrite(RHS, sizeof(scalar), ndof, f);
  fclose(f);
}

// L2 projections
template<typename Real, typename Scalar>
Scalar L2projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                           Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar L2projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                           Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (ext->fn[0]->val[i] * v->val[i]);
  return result;
}

// H1 projections
template<typename Real, typename Scalar>
Scalar H1projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                           Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar H1projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                           Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (ext->fn[0]->val[i] * v->val[i] + ext->fn[0]->dx[i] * 
                       v->dx[i] + ext->fn[0]->dy[i] * v->dy[i]);
  return result;
}

// Hcurl projections
template<typename Real, typename Scalar>
Scalar Hcurlprojection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * (u->curl[i] * conj(v->curl[i]));
    result += wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar Hcurlprojection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                              Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * (ext->fn[0]->curl[i] * conj(v->curl[i]));
    result += wt[i] * (ext->fn[0]->val0[i] * conj(v->val0[i]) + ext->fn[0]->val1[i] * conj(v->val1[i]));
  }

  return result;
}

int DiscreteProblem::assign_dofs()
{
  // sanity checks
  if (this->wf == NULL) error("this->wf = NULL in DiscreteProblem::assign_dofs().");

  // assigning dofs to each space
  if (this->spaces == NULL) error("this->spaces is NULL in DiscreteProblem::assign_dofs().");
  int ndof = 0;
  for (int i = 0; i < this->wf->neq; i++) {
    if (this->spaces[i] == NULL) error("this->spaces[%d] is NULL in assign_dofs().", i);
    int inc = this->spaces[i]->assign_dofs(ndof);
    ndof += inc;
  }

  return ndof;
}

// global orthogonal projection
void DiscreteProblem::project_global(Tuple<MeshFunction*> source, Tuple<Solution*> target, Tuple<int>proj_norms)
{
  // sanity checks
  int n = source.size();
  if (this->spaces == NULL) error("this->spaces == NULL in DiscreteProblem::project_global().");
  for (int i=0; i<n; i++) if(this->spaces[i] == NULL)
			    error("this->spaces[%d] == NULL in DiscreteProblem::project_global().", i);
  if (n != target.size())
    error("Mismatched numbers of projected functions and solutions in DiscreteProblem::project_global().");
  if (n > 10)
    error("Wrong number of projected functions in DiscreteProblem::project_global().");
  if (proj_norms != Tuple<int>()) {
    if (n != proj_norms.size())
      error("Mismatched numbers of projected functions and projection norms in DiscreteProblem::project_global().");
  }
  if (wf != NULL) {
    if (n != wf->neq)
      error("Wrong number of functions in DiscreteProblem::project_global().");
  }
  if (!have_spaces)
    error("You have to init_spaces() before using DiscreteProblem::project_global().");

  // this is needed since spaces may have their DOFs enumerated only locally
  // when they come here.
  this->assign_dofs();

  // back up original weak form
  WeakForm* wf_orig = wf;

  // define temporary projection weak form
  WeakForm wf_proj(n);
  wf = &wf_proj;
  int found[100];
  for (int i = 0; i < 100; i++) found[i] = 0;
  for (int i = 0; i < n; i++) {
    int norm;
    if (proj_norms == Tuple<int>()) norm = 1;
    else norm = proj_norms[i];
    if (norm == 0) {
      found[i] = 1;
      wf->add_matrix_form(i, i, L2projection_biform<double, scalar>, L2projection_biform<Ord, Ord>);
      wf->add_vector_form(i, L2projection_liform<double, scalar>, L2projection_liform<Ord, Ord>,
                     H2D_ANY, source[i]);
    }
    if (norm == 1) {
      found[i] = 1;
      wf->add_matrix_form(i, i, H1projection_biform<double, scalar>, H1projection_biform<Ord, Ord>);
      wf->add_vector_form(i, H1projection_liform<double, scalar>, H1projection_liform<Ord, Ord>,
                     H2D_ANY, source[i]);
    }
    if (norm == 2) {
      found[i] = 1;
      wf->add_matrix_form(i, i, Hcurlprojection_biform<double, scalar>, Hcurlprojection_biform<Ord, Ord>);
      wf->add_vector_form(i, Hcurlprojection_liform<double, scalar>, Hcurlprojection_liform<Ord, Ord>,
                     H2D_ANY, source[i]);
    }
  }
  for (int i=0; i < n; i++) {
    if (found[i] == 0) {
      printf("index of component: %d\n", i);
      error("Wrong projection norm in DiscreteProblem::project_global().");
    }
  }

  //assembling the projection matrix, Dir vector and RHS
  DiscreteProblem::assemble(this->A, this->Dir, this->RHS, false);
  // since this is a linear problem, put the Dir vector to the right-hand side:
  for (int i=0; i < this->get_num_dofs(); i++) RHS[i] += Dir[i];

  DiscreteProblem::solve(target);

  // restoring original weak form
  wf = wf_orig;
  wf_seq = -1;
}

void DiscreteProblem::project_global(Tuple<MeshFunction*> source, Tuple<Solution*> target,
                               matrix_forms_tuple_t proj_biforms, vector_forms_tuple_t proj_liforms )
{
  // sanity checks
  int n = source.size();
  if (this->spaces == NULL) error("this->spaces == NULL in DiscreteProblem::project_global().");
  for (int i=0; i<n; i++) if(this->spaces[i] == NULL)
			    error("this->spaces[%d] == NULL in DiscreteProblem::project_global().", i);
  if (n != target.size())
    error("Mismatched numbers of projected functions and solutions in DiscreteProblem::project_global().");
  if (n > 10)
    error("Wrong number of projected functions in DiscreteProblem::project_global().");

  matrix_forms_tuple_t::size_type n_biforms = proj_biforms.size();
  if (n_biforms != proj_liforms.size())
    error("Mismatched numbers of projection forms in DiscreteProblem::project_global().");
  if (n_biforms > 0) {
    if (n != n_biforms)
      error("Mismatched numbers of projected functions and projection forms in DiscreteProblem::project_global().");
  }
  else
    warn("DiscreteProblem::project_global() expected %d user defined biform(s) & liform(s); ordinary H1 projection will be performed", n);

  if (wf != NULL) {
    if (n != wf->neq)
      error("Wrong number of functions in DiscreteProblem::project_global().");
  }
  if (!have_spaces)
    error("You have to init_spaces() before using DiscreteProblem::project_global().");

  // this is needed since spaces may have their DOFs enumerated only locally
  // when they come here.
  this->assign_dofs();

  // back up original weak form
  WeakForm* wf_orig = wf;

  // define temporary projection weak form
  WeakForm wf_proj(n);
  wf = &wf_proj;
  for (int i = 0; i < n; i++) {
    if (n_biforms == 0) {
      wf->add_matrix_form(i, i, H1projection_biform<double, scalar>, H1projection_biform<Ord, Ord>);
      wf->add_vector_form(i, H1projection_liform<double, scalar>, H1projection_liform<Ord, Ord>,
                     H2D_ANY, source[i]);
    }
    else {
      wf->add_matrix_form(i, i, proj_biforms[i].first, proj_biforms[i].second);
      wf->add_vector_form(i, proj_liforms[i].first, proj_liforms[i].second,
                     H2D_ANY, source[i]);
    }
  }

  //assembling the projection matrix, Dir vector and RHS
  DiscreteProblem::assemble(this->A, this->Dir, this->RHS, false);
  // since this is a linear problem, put the Dir vector to the right-hand side:
  for (int i=0; i < this->get_num_dofs(); i++) RHS[i] += Dir[i];

  DiscreteProblem::solve(target);

  // restoring original weak form
  wf = wf_orig;
  wf_seq = -1;
}

int DiscreteProblem::get_num_dofs()
{
  // sanity checks
  if (this->wf == NULL) error("this->wf is NULL in DiscreteProblem::get_num_dofs().");
  if (this->wf->neq == 0) error("this->wf->neq is 0 in DiscreteProblem::get_num_dofs().");
  if (this->spaces == NULL) error("this->spaces[%d] is NULL in DiscreteProblem::get_num_dofs().");

  int ndof = 0;
  for (int i = 0; i < this->wf->neq; i++) {
    if (this->spaces[i] ==  NULL) error("this->spaces[%d] is NULL in DiscreteProblem::get_num_dofs().", i);
    ndof += this->get_num_dofs(i);
  }
  return ndof;
}

void DiscreteProblem::update_essential_bc_values()
{
  int n = this->wf->neq;
  for (int i=0; i<n; i++) this->spaces[i]->update_essential_bc_values();
}

// Newton's method for an arbitrary number of equations.
bool DiscreteProblem::solve_newton(Tuple<Solution*> u_prev, double newton_tol, 
                                int newton_max_iter, bool verbose, 
                                Tuple<MeshFunction*> mesh_fns) 
{
  // sanity checks
  int n = u_prev.size();
  if (n != this->wf->neq) 
    error("The number of solutions in newton_solve() must match the number of equation in the PDE system.");
  if (this->spaces == NULL) error("spaces is NULL in solve_newton().");
  for (int i=0; i < n; i++) {
    if (this->spaces[i] == NULL) error("spaces[%d] is NULL in solve_newton().", i);
  }
  int n_mesh_fns;
  if (mesh_fns == Tuple<MeshFunction*>()) n_mesh_fns = 0;
  else n_mesh_fns = mesh_fns.size();
  for (int i=0; i<n_mesh_fns; i++) {
    if (mesh_fns[i] == NULL) error("a filter is NULL in solve_newton().");
  }

  int it = 1;
  double res_l2_norm;
  do
  {
    info("---- Newton iter %d:", it); 

    // reinitialize filters
    for (int i=0; i < n_mesh_fns; i++) mesh_fns[i]->reinit();

    // assemble the Jacobian matrix and residual vector,
    // solve the system
    this->assemble(this->A, this->Dir, this->RHS, false);
    this->solve(u_prev);

    // calculate the l2-norm of residual vector
    res_l2_norm = this->get_residual_l2_norm();
    if (verbose) printf("---- Newton iter %d, ndof %d, res. l2 norm %g\n", 
                        it, this->get_num_dofs(), res_l2_norm);

    it++;
  }
  while (res_l2_norm > newton_tol && it <= newton_max_iter);

  // returning "true" if converged, otherwise returning "false"
  if (it <= newton_max_iter) return true;
  else return false;
}

