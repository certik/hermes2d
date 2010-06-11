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
#include "python_solvers.h"

void qsort_int(int* pbase, size_t total_elems); // defined in qsort.cpp

//// interface /////////////////////////////////////////////////////////////////////////////////////

void LinSystem::init(WeakForm* wf_, Solver* solver_) 
{
  this->wf = wf_;
  this->solver = solver_;
  slv_ctx = this->solver ? this->solver->new_context(false) : NULL;

  this->RHS = this->Dir = this->Vec = NULL;
  this->A = NULL;
  this->mat_sym = false;

  this->num_spaces = this->wf->neq;
  this->spaces = new Space*[this->wf->neq];
  for (int i=0; i<this->wf->neq; i++) this->spaces[i] = NULL;
  this->meshes = new Mesh*[this->wf->neq];
  for (int i=0; i<this->wf->neq; i++) this->meshes[i] = NULL;
  this->sp_seq = new int[this->wf->neq];
  this->wf_seq = -1;
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  this->pss = new PrecalcShapeset*[this->wf->neq];
  for (int i=0; i<this->wf->neq; i++) this->pss[i] = NULL;
  this->num_user_pss = 0;

  this->values_changed = true;
  this->struct_changed = true;
  this->have_spaces = false;
  this->want_dir_contrib = true;
}

// this is needed because of a constructor in NonlinSystem 
LinSystem::LinSystem() {}

LinSystem::LinSystem(WeakForm* wf, Solver* solver)
{
  this->init(wf, solver);
}

LinSystem::LinSystem(WeakForm* wf)
{
  Solver *solver = NULL;
  this->init(wf, solver);
}

LinSystem::LinSystem(WeakForm* wf, Solver* solver, Space* s)
{
  this->init(wf, solver);
  this->init_space(s);
}

LinSystem::LinSystem(WeakForm* wf, Space* s)
{
  Solver *solver = NULL;
  this->init(wf, solver);
  this->init_space(s);
}

// NOTE: there is some code duplication in the two constructors below.
// FIXME: Ivo's Tuples should be used here, then the duplication can be avoided.
LinSystem::LinSystem(WeakForm* wf, Solver* solver, int n, ...)
{
  this->init(wf, solver);
  // set spaces, meshes and pss
  if (n <= 0 || n > wf->neq) error("Bad number of spaces.");
  this->num_spaces = n;
  va_list ap;
  va_start(ap, n);
  // set spaces and meshes at the same time
  for (int i = 0; i < wf->neq; i++) {
    this->spaces[i] = (i < n) ? va_arg(ap, Space*) : this->spaces[n-1];
    this->meshes[i] = (i < n) ? this->spaces[i]->mesh : this->spaces[n-1]->mesh;
  }
  va_end(ap);
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  this->have_spaces = true;
  for (int i = 0; i < n; i++){
    if (this->pss[i] == NULL) {
      Shapeset *shapeset = this->spaces[i]->get_shapeset();
      PrecalcShapeset *p = new PrecalcShapeset(shapeset);
      if (p == NULL) error("New PrecalcShapeset could not be allocated.");
      this->pss[i] = p;
      this->num_user_pss++;
    }
  }
}

LinSystem::LinSystem(WeakForm* wf, int n, ...)
{
  Solver* solver = NULL;
  this->init(wf, solver);
  // set spaces, meshes and pss
  if (n <= 0 || n > wf->neq) error("Bad number of spaces.");
  this->num_spaces = n;
  va_list ap;
  va_start(ap, n);
  // set spaces and meshes at the same time
  for (int i = 0; i < wf->neq; i++) {
    this->spaces[i] = (i < n) ? va_arg(ap, Space*) : this->spaces[n-1];
    this->meshes[i] = (i < n) ? this->spaces[i]->mesh : this->spaces[n-1]->mesh;
  }
  va_end(ap);
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  this->have_spaces = true;
  for (int i = 0; i < n; i++){
    if (this->pss[i] == NULL) {
      Shapeset *shapeset = this->spaces[i]->get_shapeset();
      PrecalcShapeset *p = new PrecalcShapeset(shapeset);
      if (p == NULL) error("New PrecalcShapeset could not be allocated.");
      this->pss[i] = p;
      this->num_user_pss++;
    }
  }
}

LinSystem::~LinSystem()
{
  free();
  /* FIXME SOON - this is a memory leak. Uncommenting 
     this was causing a double-free segfault. */
  //free_meshes_and_spaces();
  if (this->sp_seq != NULL) delete [] this->sp_seq;
  if (this->pss != NULL) delete [] this->pss;
  if (this->solver) this->solver->free_context(this->slv_ctx);
}

void LinSystem::free_meshes_and_spaces()
{
  // free spaces, making sure that duplicated ones do not get deleted twice
  if (spaces != NULL)
  {
    for (int i = 0; i < this->num_spaces; i++) {
      for (int j = i+1; j < this->num_spaces; j++) {
        if (spaces[j] == spaces[i]) spaces[j] = NULL;
      }
    }
    for (int i = 0; i < this->num_spaces; i++) {
      if (spaces[i] != NULL) {
        delete spaces[i];
        spaces[i] = NULL;
      }
    }
    //spaces = NULL;
  }

  // free meshes, making sure that duplicated ones do not get deleted twice
  if (meshes != NULL)
  {
    for (int i = 0; i < this->num_spaces; i++) {
      for (int j = i+1; j < this->num_spaces; j++) {
        if (meshes[j] == meshes[i]) meshes[j] = NULL;
      }
    }
    for (int i = 0; i < this->num_spaces; i++) {
      if (meshes[i] != NULL) {
        delete meshes[i];
        meshes[i] = NULL;
      }
    }
    //meshes = NULL;
  }
} 

// Should not be called by the user.
void LinSystem::init_spaces(int n, ...)
{
  if (n <= 0 || n > wf->neq) error("Bad number of spaces.");
  num_spaces = n;
  va_list ap;
  va_start(ap, n);
  // set spaces and meshes at the same time
  for (int i = 0; i < wf->neq; i++) {
    spaces[i] = (i < n) ? va_arg(ap, Space*) : spaces[n-1];
    meshes[i] = (i < n) ? spaces[i]->mesh : spaces[n-1]->mesh;
  }
  va_end(ap);
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  have_spaces = true;
  for (int i = 0; i < n; i++){
    if (pss[i] == NULL) {
      Shapeset *shapeset = spaces[i]->get_shapeset();
      PrecalcShapeset *p = new PrecalcShapeset(shapeset);
      if (p == NULL) error("New PrecalcShapeset could not be allocated.");
      pss[i] = p;
      num_user_pss++;
    }
  }
}

// Deprecated.
void LinSystem::set_spaces(int n, ...)
{
  warn("Call to deprecated function set_spaces().");
  if (n <= 0 || n > wf->neq) error("Bad number of spaces.");
  num_spaces = n;
  va_list ap;
  va_start(ap, n);
  // set spaces and meshes at the same time
  for (int i = 0; i < wf->neq; i++) {
    spaces[i] = (i < n) ? va_arg(ap, Space*) : spaces[n-1];
    meshes[i] = (i < n) ? spaces[i]->mesh : spaces[n-1]->mesh;
  }
  va_end(ap);
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  have_spaces = true;
  for (int i = 0; i < n; i++){
    if (pss[i] == NULL) {
      Shapeset *shapeset = spaces[i]->get_shapeset();
      PrecalcShapeset *p = new PrecalcShapeset(shapeset);
      if (p == NULL) error("New PrecalcShapeset could not be allocated.");
      pss[i] = p;
      num_user_pss++;
   }
  }
}

// Should not be called by the user.
void LinSystem::init_space(Space* s)
{
  num_spaces = 1;
  spaces[0] = s;
  meshes[0] = s->mesh;
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  have_spaces = true;
  if (pss[0] == NULL) {
    Shapeset *shapeset = spaces[0]->get_shapeset();
    PrecalcShapeset *p = new PrecalcShapeset(shapeset);
    if (p == NULL) error("New PrecalcShapeset could not be allocated.");
    pss[0] = p;
    num_user_pss = 1;
  }
}

// Deprecated.
void LinSystem::set_space(Space* s)
{
  warn("Call to deprecated function set_space().");
  num_spaces = 1;
  spaces[0] = s;
  meshes[0] = s->mesh;
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  have_spaces = true;
  if (pss[0] == NULL) {
    Shapeset *shapeset = spaces[0]->get_shapeset();
    PrecalcShapeset *p = new PrecalcShapeset(shapeset);
    if (p == NULL) error("New PrecalcShapeset could not be allocated.");
    pss[0] = p;
    num_user_pss = 1;
  }
}

void LinSystem::set_pss(int n, ...)
{
  warn("Call to deprecated function LinSystem::set_pss().");
  if (n <= 0 || n > wf->neq) error("Bad number of pss's.");

  va_list ap;
  va_start(ap, n);
  for (int i = 0; i < n; i++)
    pss[i] = va_arg(ap, PrecalcShapeset*);
  va_end(ap);
  num_user_pss = n;

  for (int i = n; i < wf->neq; i++)
  {
    if (spaces[i]->get_shapeset() != spaces[n-1]->get_shapeset())
      error("Spaces with different shapesets must have different pss's.");
    pss[i] = new PrecalcShapeset(pss[n-1]);
  }
}

void LinSystem::set_pss(PrecalcShapeset* p)
{
  int n = 1;
  pss[0] = p;
  num_user_pss = n;

  for (int i = n; i < wf->neq; i++)
  {
    if (spaces[i]->get_shapeset() != spaces[n-1]->get_shapeset())
      error("Spaces with different shapesets must have different pss's.");
    pss[i] = new PrecalcShapeset(pss[n-1]);
  }
}

void LinSystem::copy(LinSystem* sys)
{
  error("Not implemented yet.");
}


void LinSystem::free()
{
  matrix_free();
  if (this->RHS != NULL) { ::free(this->RHS); this->RHS = NULL; }
  if (this->Dir != NULL) { ::free(this->Dir-1); this->Dir = NULL; }
  if (this->Vec != NULL) { delete [] this->Vec; this->Vec = NULL; }
  if (this->solver) this->solver->free_data(this->slv_ctx);

  this->struct_changed = this->values_changed = true;
  memset(this->sp_seq, -1, sizeof(int) * this->wf->neq);
  this->wf_seq = -1;
}

void LinSystem::matrix_free()
{
  if (this->A != NULL) { ::delete this->A; this->A = NULL; }
}

//// matrix creation ///////////////////////////////////////////////////////////////////////////////

void LinSystem::create_matrix(bool rhsonly)
{
  // sanity check
  if (!this->have_spaces)
    error("Before assemble(), you need to call init_spaces().");

  // check if we can reuse the matrix structure
  bool up_to_date = true;
  for (int i = 0; i < this->wf->neq; i++)
    if (this->spaces[i]->get_seq() != this->sp_seq[i])
      { up_to_date = false; break; }
  if (this->wf->get_seq() != this->wf_seq)
    up_to_date = false;

  // if yes, just zero the values and we're done
  if (up_to_date)
  {
    verbose("Reusing matrix sparse structure.");
    if (!rhsonly) {
#ifdef H2D_COMPLEX
      this->A = new CooMatrix(this->ndofs, true);
#else
      this->A = new CooMatrix(this->ndofs);
#endif
      memset(this->Dir, 0, sizeof(scalar) * this->ndofs);
    }
    memset(this->RHS, 0, sizeof(scalar) * this->ndofs);
    return;
  }
  else if (rhsonly)
    error("Cannot reassemble RHS only: spaces have changed.");

  // spaces have changed: create the matrix from scratch
  this->free();
  trace("Creating matrix sparse structure...");
  TimePeriod cpu_time;

  // calculate the total number of DOF
  this->ndofs = this->get_num_dofs();

  this->RHS = (scalar*) malloc(sizeof(scalar) * this->ndofs);

#ifdef H2D_COMPLEX
  this->A = new CooMatrix(this->ndofs, true);
#else
  this->A = new CooMatrix(this->ndofs);
#endif
  this->Dir = (scalar*) malloc(sizeof(scalar) * (this->ndofs + 1)) + 1;
  if (this->RHS == NULL || this->Dir == NULL) error("Out of memory. Error allocating the RHS vector.");
  memset(this->RHS, 0, sizeof(scalar) * this->ndofs);
  memset(this->Dir, 0, sizeof(scalar) * this->ndofs);

  // save space seq numbers and weakform seq number, so we can detect their changes
  for (int i = 0; i < this->wf->neq; i++)
    this->sp_seq[i] = this->spaces[i]->get_seq();
  this->wf_seq = this->wf->get_seq();

  this->struct_changed = true;
}


int LinSystem::get_matrix_size() const
{
    return this->A->get_size();
}

void LinSystem::get_matrix(int*& Ap, int*& Ai, scalar*& Ax, int& size) const
{
    /// XXX: this is a memory leak:
    CSRMatrix *m = new CSRMatrix(this->A);
    Ap = m->get_IA(); Ai = m->get_JA();
    size = m->get_size();
#ifdef H2D_COMPLEX
    Ax = m->get_A_cplx();
#else
    Ax = m->get_A();
#endif
}


//// assembly //////////////////////////////////////////////////////////////////////////////////////

void LinSystem::insert_block(scalar** mat, int* iidx, int* jidx, int ilen, int jlen)
{
    this->A->add_block(iidx, ilen, jidx, jlen, mat);
}


void LinSystem::assemble(bool rhsonly)
{
  if (!rhsonly) matrix_free();
  int k, m, n, marker;
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
  if (!ndofs) return;

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
          // we now determine this for each weakform, so we comment this out
          // for now:
        //if (!bnd[edge]) continue;
        marker = ep[edge].marker;



        // obtain the list of shape functions which are nonzero on this edge
        for (unsigned int i = 0; i < s->idx.size(); i++) {
          if (e[i] == NULL) continue;
          int j = s->idx[i];
          //nat[j] = (spaces[j]->bc_type_callback(marker) == BC_NATURAL);
          //if (nat[j])
          spaces[j]->get_edge_assembly_list(e[i], edge, &al[j]);
        }

        // assemble surface bilinear forms ///////////////////////////////////
        for (unsigned int ww = 0; ww < s->bfsurf.size(); ww++)
        {
          WeakForm::BiFormSurf* bfs = s->bfsurf[ww];
          if (isempty[bfs->i] || isempty[bfs->j]) continue;
          // now we determine whether to call the form depending on the "area"
          // and the BC "marker" parameters:
          int call_form = 0;
          if (bfs->area == H2D_ANY_EDGE) {
              // call the form on all element edges (both boundary and interior)
              call_form = 1;
          } else if (bnd[edge]) {
              if (bfs->area == H2D_ANY_BOUNDARY) {
                  // call the form on all domain boundary edges only
                  call_form = 1;
              } else if (wf->is_in_area(marker, bfs->area)) {
                  // call the form on the domain boundary edges with the
                  // correct marker only
                  call_form = 1;
              }
          }
          if (!call_form) continue;

          // old conditions:
          //if (bfs->area != H2D_ANY && !wf->is_in_area(marker, bfs->area)) continue;
          m = bfs->i;  fv = spss[m];  am = &al[m];
          n = bfs->j;  fu = pss[n];   an = &al[n];

          //if (!nat[m] || !nat[n]) continue;
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
          // now we determine whether to call the form depending on the "area"
          // and the BC "marker" parameters:
          int call_form = 0;
          if (lfs->area == H2D_ANY_EDGE) {
              // call the form on all element edges (both boundary and interior)
              call_form = 1;
          } else if (bnd[edge]) {
              if (lfs->area == H2D_ANY_BOUNDARY) {
                  // call the form on all domain boundary edges only
                  call_form = 1;
              } else if (wf->is_in_area(marker, lfs->area)) {
                  // call the form on the domain boundary edges with the
                  // correct marker only
                  call_form = 1;
              }
          }
          if (!call_form) continue;
          m = lfs->i;  fv = spss[m];  am = &al[m];

          //if (!nat[m]) continue;
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
  if (want_dir_contrib)
    for (int i = 0; i < ndofs; i++)
      RHS[i] += Dir[i];

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

  scalar res = lf->evaluate_fn(np, jwt, v, e, ext, rv->get_active_element(), fv->get_shapeset(), fv->get_active_shape());

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
  Func<scalar>** ext_fn2 = new Func<scalar>*[lf->ext.size()];
  for (int i = 0; i < lf->ext.size(); i++) {
      // This is the element, that we are currently assembling in the main
      // assembly loop
      Element *e_central = rv->get_active_element();
      // this is the "other" element, that we want to calculate the flux over
      Element *e_neigh = e_central->get_neighbor(ep->edge);
      if (e_neigh != NULL) {
          // Obtain the i-th solution from the ExtData
          MeshFunction *m = (lf->ext[i]);
          Solution *sln = dynamic_cast<Solution*>(m);

          // Obtain the corresponding edge number from the neighboring element
          Node *en = e_central->en[ep->edge];
          int edge_neigh = -1;
          for (int j=0; j < e_neigh->nvert; j++)
              if (e_neigh->en[j] == en) edge_neigh = j;
          assert(edge_neigh != -1);

          // initialize the quadrature rules
          Quad2D* quad = fv->get_quad_2d();
          int eo_neigh = quad->get_edge_points(edge_neigh);
          // obtain the solution values at the edge
          sln->set_active_element(e_neigh);
          sln->set_quad_order(eo_neigh);
          double *fn_neigh = sln->get_fn_values();
          // assign those values into ext_fn2[i]
          int np = quad->get_num_points(eo_neigh);
          int nc = fv->get_num_components();
          Func<scalar>* u = new Func<scalar>(np, nc);
          u->val = new scalar[np];
          for (int k=0; k < np; k++)
              u->val[k] = fn_neigh[k];
          ext_fn2[i] = u;
      }
  }
  ext->nf2 = lf->ext.size();
  ext->fn2 = ext_fn2;

  scalar res = lf->fn(np, jwt, v, e, ext);

  ext->free(); delete ext;
  return 0.5 * res; // Edges are parametrized from 0 to 1 while integration weights
                    // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                    // the weights.
}


//// solve /////////////////////////////////////////////////////////////////////////////////////////

bool LinSystem::solve(int n, ...)
{
  if (!this->solver) error("Cannot solve -- no matrix solver was provided.");
  TimePeriod cpu_time;

  // solve the system
  if (this->Vec != NULL) delete [] this->Vec;
  this->Vec = new scalar[this->ndofs];
  memcpy(this->Vec, this->RHS, sizeof(scalar) * this->ndofs);
  solve_linear_system_scipy_umfpack(this->A, this->Vec);
  report_time("LinSystem solved in %g s", cpu_time.tick().last());

  // initialize the Solution classes
  va_list ap;
  va_start(ap, n);
  if (n > this->wf->neq) n = this->wf->neq;
  for (int i = 0; i < n; i++)
  {
    Solution* sln = va_arg(ap, Solution*);
    sln->set_fe_solution(this->spaces[i], this->pss[i], this->Vec);
  }
  va_end(ap);
  report_time("Exported solution in %g s", cpu_time.tick().last());
  return true;
}

bool LinSystem::solve(Solution* sln)
{
  if (!this->solver) error("Cannot solve -- no solver was provided.");
  TimePeriod cpu_time;

  // solve the system
  if (this->Vec != NULL) delete [] this->Vec;
  //this->Vec = (scalar*) malloc(this->ndofs * sizeof(scalar));
  this->Vec = new scalar[this->ndofs];
  if (this->Vec == NULL) error("Malloc problem in LinSystem::solve()");

  // copy the content of RHS into Vec
  memcpy(this->Vec, this->RHS, sizeof(scalar) * this->ndofs);

  // call UMFpack in Scipy
  solve_linear_system_scipy_umfpack(this->A, this->Vec);
  report_time("LinSystem solved in %g s", cpu_time.tick().last());

  // initialize the Solution class
  if (this->spaces[0] == NULL) error("No spaces in LinSystem::solve().");
  if (this->pss[0] == NULL) error("No precalc shapeset in LinSystem::solve().");
  sln->set_fe_solution(this->spaces[0], this->pss[0], this->Vec);
  report_time("Exported solution in %g s", cpu_time.tick().last());

  return true;
}


//// matrix and solution output /////////////////////////////////////////////////////////////////////////////////

void LinSystem::get_solution_vector(std::vector<scalar>& sln_vector_out) const {
  assert_msg(ndofs > 0, "Number of DOFs is not greater than zero");
  std::vector<scalar> temp(ndofs);
  sln_vector_out.swap(temp);
  for(int i = 0; i < ndofs; i++)
    sln_vector_out[i] = Vec[i];
}

void LinSystem::save_matrix_matlab(const char* filename, const char* varname)
{
}


void LinSystem::save_rhs_matlab(const char* filename, const char* varname)
{
  if (RHS == NULL) error("RHS has not been assembled yet.");
  FILE* f = fopen(filename, "w");
  if (f == NULL) error("Could not open file %s for writing.", filename);
  verbose("Saving RHS vector in MATLAB format...");
  fprintf(f, "%% Size: %dx1\n%s = [\n", ndofs, varname);
  for (int i = 0; i < ndofs; i++)
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
  if (RHS == NULL) error("RHS has not been assembled yet.");
  FILE* f = fopen(filename, "wb");
  if (f == NULL) error("Could not open file %s for writing.", filename);
  verbose("Saving RHS vector in binary format...");
  hermes2d_fwrite("H2DR\001\000\000\000", 1, 8, f);
  int ssize = sizeof(scalar);
  hermes2d_fwrite(&ssize, sizeof(int), 1, f);
  hermes2d_fwrite(&ndofs, sizeof(int), 1, f);
  hermes2d_fwrite(RHS, sizeof(scalar), ndofs, f);
  fclose(f);
}

void LinSystem::set_vec_zero()
{
  if (this->Vec != NULL) delete [] this->Vec;
  this->Vec = new scalar[this->ndofs];
  for(int i=0; i<this->ndofs; i++) this->Vec[i] = 0;
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

void LinSystem::project_global_n(int proj_norm, int n, ...)
{
  if (!have_spaces)
    error("You have to init_spaces() before using project_global_n().");
  if (n != wf->neq || n > 10)
    error("Wrong number of functions in project_global_n().");

  int i;
  MeshFunction* fn[10];
  Solution* result[10];

  va_list ap;
  va_start(ap, n);
  for (i = 0; i < n; i++)
    fn[i] = va_arg(ap, MeshFunction*);
  for (i = 0; i < n; i++)
    result[i] = va_arg(ap, Solution*);
  va_end(ap);

  WeakForm* wf_orig = wf;

  WeakForm wf_proj(n);
  wf = &wf_proj;
  int found = 0;
  if (proj_norm == 0) {
    found = 1;
    for (i = 0; i < n; i++)
    {
      wf->add_biform(i, i, L2projection_biform<double, scalar>, L2projection_biform<Ord, Ord>);
      wf->add_liform(i, L2projection_liform<double, scalar>, L2projection_liform<Ord, Ord>, H2D_ANY, 1, fn[i]);
    }
  }
  if (proj_norm == 1) {
    found = 1;
    for (i = 0; i < n; i++)
    {
      wf->add_biform(i, i, H1projection_biform<double, scalar>, H1projection_biform<Ord, Ord>);
      wf->add_liform(i, H1projection_liform<double, scalar>, H1projection_liform<Ord, Ord>, H2D_ANY, 1, fn[i]);
    }
  }
  if (proj_norm == 2) {
    found = 1;
    for (i = 0; i < n; i++)
    {
      wf->add_biform(i, i, Hcurlprojection_biform<double, scalar>, Hcurlprojection_biform<Ord, Ord>);
      wf->add_liform(i, Hcurlprojection_liform<double, scalar>, Hcurlprojection_liform<Ord, Ord>, H2D_ANY, 1, fn[i]);
    }
  }
  if (!found) error("Wrong projection type in project_global_n().");

  want_dir_contrib = true;
  LinSystem::assemble();
  LinSystem::solve(n, result[0], result[1], result[2], result[3], result[4],
                      result[5], result[6], result[7], result[8], result[9]);
  want_dir_contrib = false;

  wf = wf_orig;
  wf_seq = -1;
}

int LinSystem::get_num_dofs()
{
  int ndof = 0;
  for (int i = 0; i < this->num_spaces; i++) {
    ndof += this->spaces[i]->get_num_dofs();
  }
  if (ndof == 0)
    error("Zero matrix size while creating matrix sparse structure.");
  return ndof;
}
