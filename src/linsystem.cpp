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
#include "linsystem.h"
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

void LinSystem::init_lin(WeakForm* wf_, CommonSolver* solver_)
{
  if (wf_ == NULL) error("LinSystem: a weak form must be given.");
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
  this->want_dir_contrib = true;
}

// this is needed because of a constructor in NonlinSystem
LinSystem::LinSystem() {}

LinSystem::LinSystem(WeakForm* wf_, CommonSolver* solver_)
{
  this->init_lin(wf_, solver_);
}

LinSystem::LinSystem(WeakForm* wf_)
{
  CommonSolver *solver_ = NULL;
  this->init_lin(wf_, solver_);
}

LinSystem::LinSystem(WeakForm* wf_, CommonSolver* solver_, Tuple<Space*> sp)
{
  int n = sp.size();
  if (wf_ == NULL) warn("Weak form is NULL.");
  if (wf_ != NULL) {
    if (n != wf_->neq)
      error("Number of spaces does not match number of equations in LinSystem::LinSystem().");
  }
  this->init_lin(wf_, solver_);
  this->init_spaces(sp);
  this->alloc_and_zero_vectors();
}

LinSystem::LinSystem(WeakForm* wf_, Tuple<Space*> sp)
{
  CommonSolver* solver_ = NULL;
  this->init_lin(wf_, solver_);
  this->init_spaces(sp);
  this->alloc_and_zero_vectors();
}

LinSystem::LinSystem(WeakForm* wf_, CommonSolver* solver_, Space* s_)
{
  if (wf_ == NULL) warn("Weak form is NULL.");
  if (wf_ != NULL) {
    if (wf_->neq != 1)
      error("Number of spaces does not match number of equations in LinSystem::LinSystem().");
  }
  this->init_lin(wf_, solver_);
  this->init_space(s_);
  this->alloc_and_zero_vectors();
}

LinSystem::LinSystem(WeakForm* wf_, Space* s_)
{
  CommonSolver *solver_ = NULL;
  this->init_lin(wf_, solver_);
  this->init_space(s_);
  this->alloc_and_zero_vectors();
}

LinSystem::LinSystem(WeakForm* wf_, CommonSolver* solver_, Space* space1_, Space* space2_)
{
  int n = 2;
  if (wf_ == NULL) warn("Weak form is NULL.");
  if (wf_ != NULL) {
    if (n != wf_->neq)
      error("Number of spaces does not match number of equations in LinSystem::LinSystem().");
  }
  this->init_lin(wf_, solver_);
  this->init_spaces(Tuple<Space*>(space1_, space2_));
  this->alloc_and_zero_vectors();
}

LinSystem::~LinSystem()
{
  /* FIXME - this should be uncommented but then it gives double-free
             segfaults in adaptive examples.
  free();
  if (this->sp_seq != NULL) delete [] this->sp_seq;
  if (this->pss != NULL) delete [] this->pss;
  if (this->solver != NULL) this->solver->free_context(this->slv_ctx);
  */
  free_vectors();
}

void LinSystem::free_spaces()
{
  // free spaces, making sure that duplicated ones do not get deleted twice
  if (this->spaces != NULL)
  {
    // this loop skipped if there is only one space
    for (int i = 0; i < this->wf->neq; i++) {
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
void LinSystem::init_spaces(Tuple<Space*> sp)
{
  int n = sp.size();
  if (n != this->wf->neq)
    error("Number of spaces does not match number of equations in LinSystem::init_spaces().");

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
    if (shapeset == NULL) error("Internal in LinSystem::init_spaces().");
    PrecalcShapeset *p = new PrecalcShapeset(shapeset);
    if (p == NULL) error("New PrecalcShapeset could not be allocated in LinSystem::init_spaces().");
    this-> pss[i] = p;
    this->num_user_pss++;
  }
}

// Should not be called by the user.
void LinSystem::init_space(Space* s)
{
  if (this->wf->neq != 1)
    error("Do not call init_space() for PDE systems, call init_spaces() instead.");
  this->init_spaces(Tuple<Space*>(s));
}

// Obsolete. Should be removed after FeProblem is removed.
void LinSystem::set_spaces(Tuple<Space*>spaces)
{
  this->init_spaces(spaces);
}

void LinSystem::set_pss(Tuple<PrecalcShapeset*> pss)
{
  warn("Call to deprecated function LinSystem::set_pss().");
  int n = pss.size();
  if (n != this->wf->neq)
    error("The number of precalculated shapesets must match the number of equations.");

  for (int i = 0; i < n; i++) this->pss[i] = pss[i];
  num_user_pss = n;
}

void LinSystem::set_pss(PrecalcShapeset* pss)
{
  this->set_pss(Tuple<PrecalcShapeset*>(pss));
}

void LinSystem::copy(LinSystem* sys)
{
  error("Not implemented yet.");
}

void LinSystem::free_vectors()
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

void LinSystem::alloc_and_zero_vectors()
{
  int ndof = this->get_num_dofs();
  //printf("debug: allocating vectors Vec, RHS, Dir for ndof   %d\n", ndof);

  if (this->RHS != NULL || this->Dir != NULL || this->Vec != NULL)
    error("All vectors must be NULL in alloc_and_zero_vectors() to prevent loss of information.");

  this->Vec = new scalar[ndof];
  if (Vec == NULL) error("Not enough memory LinSystem::alloc_and_zero_vectors().");
  memset(this->Vec, 0, ndof*sizeof(scalar));
  this->Vec_length = ndof;

  this->RHS = new scalar[ndof];
  if (RHS == NULL) error("Not enough memory in LinSystem::alloc_and_zero_vectors().");
  memset(this->RHS, 0, ndof*sizeof(scalar));
  this->RHS_length = ndof;

  this->Dir = new scalar[ndof];
  if (Dir == NULL) error("Not enough memory in LinSystem::alloc_and_zero_vectors().");
  memset(this->Dir, 0, ndof*sizeof(scalar));
  this->Dir_length = ndof;
}

void LinSystem::realloc_and_zero_vectors()
{
  int ndof = this->get_num_dofs();
  //printf("debug: reallocating vectors Vec, RHS, Dir length   %d -> %d\n", this->Vec_length, ndof);

  this->Vec = (scalar*)realloc(this->Vec, ndof*sizeof(scalar));
  if (this->Vec == NULL) error("Not enough memory LinSystem::realloc_and_zero_vectors().");
  memset(this->Vec, 0, ndof*sizeof(scalar));
  this->Vec_length = ndof;

  this->RHS = (scalar*)realloc(this->RHS, ndof*sizeof(scalar));
  if (this->RHS == NULL) error("Not enough memory LinSystem::realloc_and_zero_vectors().");
  memset(this->RHS, 0, ndof*sizeof(scalar));
  this->RHS_length = ndof;

  this->Dir = (scalar*)realloc(this->Dir, ndof*sizeof(scalar));
  if (this->Dir == NULL) error("Not enough memory LinSystem::realloc_and_zero_vectors().");
  memset(this->Dir, 0, ndof*sizeof(scalar));
  this->Dir_length = ndof;
}

void LinSystem::free()
{
  free_matrix();
  free_vectors();
  free_spaces();
  delete solver_default;

  this->struct_changed = this->values_changed = true;
  memset(this->sp_seq, -1, sizeof(int) * this->wf->neq);
  this->wf_seq = -1;
}

void LinSystem::free_matrix()
{
  if (this->A != NULL) { ::delete this->A; this->A = NULL; }
}

//// matrix creation ///////////////////////////////////////////////////////////////////////////////

void LinSystem::create_matrix(bool rhsonly)
{
  // sanity checks
  if (this->wf == NULL) error("this->wf is NULL in LinSystem::get_num_dofs().");
  if (this->wf->neq == 0) error("this->wf->neq is 0 in LinSystem::get_num_dofs().");
  if (this->spaces == NULL) error("this->spaces[%d] is NULL in LinSystem::get_num_dofs().");

  // check if we can reuse the matrix structure
  bool up_to_date = true;
  for (int i = 0; i < this->wf->neq; i++) {
    if (this->spaces[i] ==  NULL) error("this->spaces[%d] is NULL in LinSystem::get_num_dofs().", i);
    if (this->spaces[i]->get_seq() != this->sp_seq[i]) {
      up_to_date = false;
      break;
    }
  }
  if (this->wf->get_seq() != this->wf_seq) up_to_date = false;

  // calculate the number of DOF
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof = 0 in LinSystem::create_matrix().");

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


int LinSystem::get_matrix_size()
{
    return this->A->get_size();
}

void LinSystem::get_matrix(int*& Ap, int*& Ai, scalar*& Ax, int& size)
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

void LinSystem::insert_block(scalar** mat, int* iidx, int* jidx, int ilen, int jlen)
{
    this->A->add_block(iidx, ilen, jidx, jlen, mat);
}

void LinSystem::assemble(bool rhsonly)
{
  // sanity checks
  if (this->have_spaces == false)
    error("Before assemble(), you need to initialize spaces.");
  if (this->wf == NULL) error("this->wf = NULL in LinSystem::assemble().");
  if (this->spaces == NULL) error("this->spaces = NULL in LinSystem::assemble().");
  int n = this->wf->neq;
  for (int i=0; i<n; i++) if (this->spaces[i] == NULL)
			    error("this->spaces[%d] is NULL in LinSystem::assemble().", i);

  // enumerate DOF to get new length of the vectors Vec, RHS and Dir,
  // and realloc these vectors if needed
  this->assign_dofs();
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof = 0 in LinSystem::assemble().");
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
  // In such a case, the bilinear forms are assembled over one mesh, and only the rhs
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
      //// assemble volume bilinear forms //////////////////////////////////////
      for (unsigned int ww = 0; ww < s->bfvol.size(); ww++)
      {
        WeakForm::BiFormVol* bfv = s->bfvol[ww];
        if (isempty[bfv->i] || isempty[bfv->j]) continue;
        if (bfv->area != H2D_ANY && !wf->is_in_area(marker, bfv->area)) continue;
        m = bfv->i;  fv = spss[m];  am = &al[m];
        n = bfv->j;  fu = pss[n];   an = &al[n];
        bool tra = (m != n) && (bfv->sym != 0);
        bool sym = (m == n) && (bfv->sym == 1);

        // assemble the local stiffness matrix for the form bfv
        scalar bi, **mat = get_matrix_buffer(std::max(am->cnt, an->cnt));
        for (int i = 0; i < am->cnt; i++)
        {
          if (!tra && (k = am->dof[i]) < 0) continue;
          fv->set_active_shape(am->idx[i]);

          if (!sym) // unsymmetric block
          {
            for (int j = 0; j < an->cnt; j++) {
              fu->set_active_shape(an->idx[j]);
              bi = eval_form(bfv, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
              if (an->dof[j] < 0) Dir[k] -= bi; else mat[i][j] = bi;
            }
          }
          else // symmetric block
          {
            for (int j = 0; j < an->cnt; j++) {
              if (j < i && an->dof[j] >= 0) continue;
              fu->set_active_shape(an->idx[j]);
              bi = eval_form(bfv, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
              if (an->dof[j] < 0) Dir[k] -= bi; else mat[i][j] = mat[j][i] = bi;
            }
          }
        }

        // insert the local stiffness matrix into the global one
        insert_block(mat, am->dof, an->dof, am->cnt, an->cnt);

        // insert also the off-diagonal (anti-)symmetric block, if required
        if (tra)
        {
          if (bfv->sym < 0) chsgn(mat, am->cnt, an->cnt);
          transpose(mat, am->cnt, an->cnt);
          insert_block(mat, an->dof, am->dof, an->cnt, am->cnt);

          // we also need to take care of the RHS...
          for (int j = 0; j < am->cnt; j++)
            if (am->dof[j] < 0)
              for (int i = 0; i < an->cnt; i++)
                if (an->dof[i] >= 0)
                  Dir[an->dof[i]] -= mat[i][j];
        }
      }

      //// assemble volume linear forms ////////////////////////////////////////
      for (unsigned int ww = 0; ww < s->lfvol.size(); ww++)
      {
        WeakForm::LiFormVol* lfv = s->lfvol[ww];
        if (isempty[lfv->i]) continue;
        if (lfv->area != H2D_ANY && !wf->is_in_area(marker, lfv->area)) continue;
        m = lfv->i;  fv = spss[m];  am = &al[m];

        for (int i = 0; i < am->cnt; i++)
        {
          if (am->dof[i] < 0) continue;
          fv->set_active_shape(am->idx[i]);
          RHS[am->dof[i]] += eval_form(lfv, fv, &refmap[m]) * am->coef[i];
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

        // assemble surface bilinear forms ///////////////////////////////////
        for (unsigned int ww = 0; ww < s->bfsurf.size(); ww++)
        {
          WeakForm::BiFormSurf* bfs = s->bfsurf[ww];
          if (isempty[bfs->i] || isempty[bfs->j]) continue;
          if (bfs->area != H2D_ANY && !wf->is_in_area(marker, bfs->area)) continue;
          m = bfs->i;  fv = spss[m];  am = &al[m];
          n = bfs->j;  fu = pss[n];   an = &al[n];

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
              bi = eval_form(bfs, fu, fv, &refmap[n], &refmap[m], &(ep[edge])) * an->coef[j] * am->coef[i];
              if (an->dof[j] >= 0) mat[i][j] = bi; else Dir[k] -= bi;
            }
          }
          insert_block(mat, am->dof, an->dof, am->cnt, an->cnt);
        }

        // assemble surface linear forms /////////////////////////////////////
        for (unsigned int ww = 0; ww < s->lfsurf.size(); ww++)
        {
          WeakForm::LiFormSurf* lfs = s->lfsurf[ww];
          if (isempty[lfs->i]) continue;
          if (lfs->area != H2D_ANY && !wf->is_in_area(marker, lfs->area)) continue;
          m = lfs->i;  fv = spss[m];  am = &al[m];

          if (!nat[m]) continue;
          ep[edge].base = trav.get_base();
          ep[edge].space_v = spaces[m];

          for (int i = 0; i < am->cnt; i++)
          {
            if (am->dof[i] < 0) continue;
            fv->set_active_shape(am->idx[i]);
            RHS[am->dof[i]] += eval_form(lfs, fv, &refmap[m], &(ep[edge])) * am->coef[i];
          }
        }
      }
      delete_cache();
    }
    trav.finish();
  }

  // add to RHS the dirichlet contributions
  if (want_dir_contrib) {
    for (int i = 0; i < ndof; i++) {
      this->RHS[i] += this->Dir[i];
    }
  }

  verbose("Stiffness matrix assembled (stages: %d)", stages.size());
  report_time("Stiffness matrix assembled in %g s", cpu_time.tick().last());
  for (int i = 0; i < wf->neq; i++) delete spss[i];
  delete [] buffer;

  if (!rhsonly) values_changed = true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Initialize integration order for external functions
ExtData<Ord>* LinSystem::init_ext_fns_ord(std::vector<MeshFunction *> &ext)
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
ExtData<scalar>* LinSystem::init_ext_fns(std::vector<MeshFunction *> &ext, RefMap *rm, const int order)
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
Func<double>* LinSystem::get_fn(PrecalcShapeset *fu, RefMap *rm, const int order)
{
  Key key(256 - fu->get_active_shape(), order, fu->get_transform(), fu->get_shapeset()->get_id());
  if (cache_fn[key] == NULL)
    cache_fn[key] = init_fn(fu, rm, order);

  return cache_fn[key];
}

// Caching transformed values
void LinSystem::init_cache()
{
  for (int i = 0; i < g_max_quad + 1 + 4; i++)
  {
    cache_e[i] = NULL;
    cache_jwt[i] = NULL;
  }
}

void LinSystem::delete_cache()
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

////////////////////////////////////////////////////////////////////////////////////////////////////

// Actual evaluation of volume bilinear form (calculates integral)
scalar LinSystem::eval_form(WeakForm::BiFormVol *bf, PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv)
{
  // determine the integration order
  int inc = (fu->get_num_components() == 2) ? 1 : 0;
  Func<Ord>* ou = init_fn_ord(fu->get_fn_order() + inc);
  Func<Ord>* ov = init_fn_ord(fv->get_fn_order() + inc);
  ExtData<Ord>* fake_ext = init_ext_fns_ord(bf->ext);

  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  Ord o = bf->ord(1, &fake_wt, ou, ov, fake_e, fake_ext);
  int order = ru->get_inv_ref_order();
  order += o.get_order();
  limit_order(order);

  ou->free_ord(); delete ou;
  ov->free_ord(); delete ov;
  delete fake_e;
  fake_ext->free_ord(); delete fake_ext;

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
  Func<double>* u = get_fn(fu, ru, order);
  Func<double>* v = get_fn(fv, rv, order);
  ExtData<scalar>* ext = init_ext_fns(bf->ext, rv, order);

  scalar res = bf->fn(np, jwt, u, v, e, ext);
  ext->free(); delete ext;
  return res;
}


// Actual evaluation of volume linear form (calculates integral)
scalar LinSystem::eval_form(WeakForm::LiFormVol *lf, PrecalcShapeset *fv, RefMap *rv)
{
  // determine the integration order
  int inc = (fv->get_num_components() == 2) ? 1 : 0;
  Func<Ord>* ov = init_fn_ord(fv->get_fn_order() + inc);
  ExtData<Ord>* fake_ext = init_ext_fns_ord(lf->ext);

  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  Ord o = lf->evaluate_ord(1, &fake_wt, ov, fake_e, fake_ext, rv->get_active_element(), fv->get_shapeset(), fv->get_active_shape());
  int order = rv->get_inv_ref_order();
  order += o.get_order();
  limit_order(order);

  ov->free_ord(); delete ov;
  delete fake_e;
  fake_ext->free_ord(); delete fake_ext;

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
  Func<double>* v = get_fn(fv, rv, order);
  ExtData<scalar>* ext = init_ext_fns(lf->ext, rv, order);

  scalar res = lf->evaluate_fn(np, jwt, v, e, ext, rv->get_active_element(),
               fv->get_shapeset(), fv->get_active_shape());

  ext->free(); delete ext;
  return res;

}


// Actual evaluation of surface bilinear form (calculates integral)
scalar LinSystem::eval_form(WeakForm::BiFormSurf *bf, PrecalcShapeset *fu,
                            PrecalcShapeset *fv, RefMap *ru, RefMap *rv, EdgePos* ep)
{
  // eval the form
  Quad2D* quad = fu->get_quad_2d();
  assert(ep->edge < 5);
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
  Func<double>* u = get_fn(fu, ru, eo);
  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(bf->ext, rv, eo);

  scalar res = bf->fn(np, jwt, u, v, e, ext);

  ext->free(); delete ext;
  return 0.5 * res; // Edges are parameterized from 0 to 1 while integration weights
                    // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                    // the weights.
}


// Actual evaluation of surface linear form (calculates integral)
scalar LinSystem::eval_form(WeakForm::LiFormSurf *lf, PrecalcShapeset *fv, RefMap *rv, EdgePos* ep)
{
  // eval the form
  Quad2D* quad = fv->get_quad_2d();
  assert(ep->edge < 5);
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
  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(lf->ext, rv, eo);

  scalar res = lf->fn(np, jwt, v, e, ext);

  ext->free(); delete ext;
  return 0.5 * res; // Edges are parametrized from 0 to 1 while integration weights
                    // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                    // the weights.
}


//// solve /////////////////////////////////////////////////////////////////////////////////////////

bool LinSystem::solve(Tuple<Solution*> sln)
{
  int n = sln.size();

  // if the number of solutions does not match the number of equations, throw error
  if (n != this->wf->neq)
    error("Number of solutions does not match the number of equations in LinSystem::solve().");

  // if Vec is not initialized, throw error
  if (this->Vec == NULL) error("Vec is NULL in LinSystem::solve().");

  // check vector size
  int ndof = this->get_num_dofs();
  if (this->Vec_length != ndof || this->RHS_length != ndof || this->Dir_length != ndof)
    error("Length of vectors Vec, RHS or Dir does not match this->ndof in LinSystem::solve().");

  // check matrix size
  if (ndof == 0) error("ndof = 0 in LinSystem::solve().");
  if (ndof != this->A->get_size())
    error("Matrix size does not match vector length in LinSystem:solve().");

  // time measurement
  TimePeriod cpu_time;

  // solve the system
  memcpy(this->Vec, this->RHS, sizeof(scalar) * ndof);
  this->solver->solve(this->A, this->Vec);
  report_time("LinSystem solved in %g s", cpu_time.tick().last());

  // copy solution coefficient vectors into Solutions
  for (int i = 0; i < n; i++)
  {
    sln[i]->set_fe_solution(this->spaces[i], this->pss[i], this->Vec);
  }

  report_time("Exported solution in %g s", cpu_time.tick().last());
  return true;
}

// single equation case
bool LinSystem::solve(Solution* sln)
{
  bool flag;
  flag = this->solve(Tuple<Solution*>(sln));
  return flag;
}

// two equations case
bool LinSystem::solve(Solution* sln1, Solution *sln2)
{
  bool flag;
  flag = this->solve(Tuple<Solution*>(sln1, sln2));
  return flag;
}

// three equations case
bool LinSystem::solve(Solution* sln1, Solution *sln2, Solution* sln3)
{
  bool flag;
  flag = this->solve(Tuple<Solution*>(sln1, sln2, sln3));
  return flag;
}


//// matrix and solution output /////////////////////////////////////////////////////////////////////////////////

void LinSystem::get_solution_vector(std::vector<scalar>& sln_vector_out) {
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof = 0 in LinSystem::get_solution_vector().");

  std::vector<scalar> temp(ndof);
  sln_vector_out.swap(temp);
  for(int i = 0; i < ndof; i++)
    sln_vector_out[i] = Vec[i];
}

void LinSystem::save_matrix_matlab(const char* filename, const char* varname)
{
  warn("Saving matrix in Matlab format not implemented yet.");
}

void LinSystem::save_rhs_matlab(const char* filename, const char* varname)
{
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof = 0 in LinSystem::save_rhs_matlab().");
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


void LinSystem::save_matrix_bin(const char* filename)
{
}

void LinSystem::save_rhs_bin(const char* filename)
{
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof = 0 in LinSystem::save_rhs_bin().");
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
Scalar L2projection_biform(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar L2projection_liform(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (ext->fn[0]->val[i] * v->val[i]);
  return result;
}

// H1 projections
template<typename Real, typename Scalar>
Scalar H1projection_biform(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar H1projection_liform(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (ext->fn[0]->val[i] * v->val[i] + ext->fn[0]->dx[i] * v->dx[i] + ext->fn[0]->dy[i] * v->dy[i]);
  return result;
}

// Hcurl projections
template<typename Real, typename Scalar>
Scalar Hcurlprojection_biform(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * (u->curl[i] * conj(v->curl[i]));
    result += wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar Hcurlprojection_liform(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * (ext->fn[0]->curl[i] * conj(v->curl[i]));
    result += wt[i] * (ext->fn[0]->val0[i] * conj(v->val0[i]) + ext->fn[0]->val1[i] * conj(v->val1[i]));
  }

  return result;
}

int LinSystem::assign_dofs()
{
  // sanity checks
  if (this->wf == NULL) error("this->wf = NULL in LinSystem::assign_dofs().");

  // assigning dofs to each space
  if (this->spaces == NULL) error("this->spaces is NULL in LinSystem::assign_dofs().");
  int ndof = 0;
  for (int i = 0; i < this->wf->neq; i++) {
    if (this->spaces[i] == NULL) error("this->spaces[%d] is NULL in assign_dofs().", i);
    int inc = this->spaces[i]->assign_dofs(ndof);
    ndof += inc;
  }

  return ndof;
}

// global orthogonal projection
void LinSystem::project_global(Tuple<MeshFunction*> source, Tuple<Solution*> target, Tuple<int>proj_norms)
{
  // sanity checks
  int n = source.size();
  if (this->spaces == NULL) error("this->spaces == NULL in LinSystem::project_global().");
  for (int i=0; i<n; i++) if(this->spaces[i] == NULL)
			    error("this->spaces[%d] == NULL in LinSystem::project_global().", i);
  if (n != target.size())
    error("Mismatched numbers of projected functions and solutions in LinSystem::project_global().");
  if (n > 10)
    error("Wrong number of projected functions in LinSystem::project_global().");
  if (proj_norms != Tuple<int>()) {
    if (n != proj_norms.size())
      error("Mismatched numbers of projected functions and projection norms in LinSystem::project_global().");
  }
  if (wf != NULL) {
    if (n != wf->neq)
      error("Wrong number of functions in LinSystem::project_global().");
  }
  if (!have_spaces)
    error("You have to init_spaces() before using LinSystem::project_global().");

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
      error("Wrong projection norm in LinSystem::project_global().");
    }
  }

  want_dir_contrib = true;
  LinSystem::assemble();
  LinSystem::solve(target);
  want_dir_contrib = false;

  // restoring original weak form
  wf = wf_orig;
  wf_seq = -1;
}

void LinSystem::project_global( Tuple<MeshFunction*> source, Tuple<Solution*> target,
                                biforms_tuple_t proj_biforms, liforms_tuple_t proj_liforms )
{
  // sanity checks
  int n = source.size();
  if (this->spaces == NULL) error("this->spaces == NULL in LinSystem::project_global().");
  for (int i=0; i<n; i++) if(this->spaces[i] == NULL)
			    error("this->spaces[%d] == NULL in LinSystem::project_global().", i);
  if (n != target.size())
    error("Mismatched numbers of projected functions and solutions in LinSystem::project_global().");
  if (n > 10)
    error("Wrong number of projected functions in LinSystem::project_global().");

  biforms_tuple_t::size_type n_biforms = proj_biforms.size();
  if (n_biforms != proj_liforms.size())
    error("Mismatched numbers of projection forms in LinSystem::project_global().");
  if (n_biforms > 0) {
    if (n != n_biforms)
      error("Mismatched numbers of projected functions and projection forms in LinSystem::project_global().");
  }
  else
    warn("LinSystem::project_global() expected %d user defined biform(s) & liform(s); ordinary H1 projection will be performed", n);

  if (wf != NULL) {
    if (n != wf->neq)
      error("Wrong number of functions in LinSystem::project_global().");
  }
  if (!have_spaces)
    error("You have to init_spaces() before using LinSystem::project_global().");

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

  want_dir_contrib = true;
  LinSystem::assemble();
  LinSystem::solve(target);
  want_dir_contrib = false;

  // restoring original weak form
  wf = wf_orig;
  wf_seq = -1;
}

int LinSystem::get_num_dofs()
{
  // sanity checks
  if (this->wf == NULL) error("this->wf is NULL in LinSystem::get_num_dofs().");
  if (this->wf->neq == 0) error("this->wf->neq is 0 in LinSystem::get_num_dofs().");
  if (this->spaces == NULL) error("this->spaces[%d] is NULL in LinSystem::get_num_dofs().");

  int ndof = 0;
  for (int i = 0; i < this->wf->neq; i++) {
    if (this->spaces[i] ==  NULL) error("this->spaces[%d] is NULL in LinSystem::get_num_dofs().", i);
    ndof += this->get_num_dofs(i);
  }
  return ndof;
}

void LinSystem::update_essential_bc_values()
{
  int n = this->wf->neq;
  for (int i=0; i<n; i++) this->spaces[i]->update_essential_bc_values();
}
