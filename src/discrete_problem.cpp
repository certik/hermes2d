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
}

DiscreteProblem::DiscreteProblem(WeakForm* wf_, Space* s_)
{
  CommonSolver *solver_ = NULL;
  this->init(wf_, solver_);
  this->init_space(s_);
}

DiscreteProblem::~DiscreteProblem()
{
  // FIXME: more things should be deleted here
  // to avoid memory leaks.
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

void DiscreteProblem::free()
{
  free_spaces();

  this->struct_changed = this->values_changed = true;
  memset(this->sp_seq, -1, sizeof(int) * this->wf->neq);
  this->wf_seq = -1;
}

//// assembly //////////////////////////////////////////////////////////////////////////////////////

void DiscreteProblem::insert_block(Matrix *mat_ext, scalar** mat, int* iidx, int* jidx, int ilen, int jlen)
{
    mat_ext->add_block(iidx, ilen, jidx, jlen, mat);
}

void DiscreteProblem::assemble(Matrix* mat_ext, Vector* dir_ext, Vector* rhs_ext, bool rhsonly)
{
  // sanity checks
  if (this->have_spaces == false)
    error("Before assemble(), you need to initialize spaces.");
  if (this->wf == NULL) error("this->wf = NULL in DiscreteProblem::assemble().");
  if (this->spaces == NULL) error("this->spaces = NULL in DiscreteProblem::assemble().");
  int n = this->wf->neq;
  for (int i = 0; i < n; i++) {
    if (this->spaces[i] == NULL) error("this->spaces[%d] is NULL in DiscreteProblem::assemble().", i);
  }
  if (rhs_ext == NULL) error("rhs_ext == NULL in DiscreteProblem::assemble().");
  if (rhsonly == false) {
    if (mat_ext == NULL) error("mat_ext == NULL in DiscreteProblem::assemble().");
    if (mat_ext->get_size() != rhs_ext->get_size()) {
	printf("mat_ext matrix size = %d\n", mat_ext->get_size());
	printf("rhs_ext vector size = %d\n", rhs_ext->get_size());
        error("Mismatched mat_ext and rhs_ext vector sizes in DiscreteProblem::assemble().");
    }
  }

  // Assign dof in all spaces. 
  // Realloc mat_ext, dir_ext and rhs_ext if ndof changed, 
  // and clear dir_ext and rhs_ext. 
  //Do not touch the matrix if rhsonly == true. 
  int ndof = this->assign_dofs();
  //printf("ndof = %d\n", ndof);
  if (ndof == 0) error("ndof = 0 in DiscreteProblem::assemble().");
  if (rhsonly == false) {
    if (mat_ext->get_size() != ndof) {
      mat_ext->init();
    }
  }
  if (dir_ext != NULL) {
    if (dir_ext->get_size() != ndof) {
      dir_ext->free_data();
      dir_ext->init(ndof);
    }
    else dir_ext->set_zero();
  }
  if (rhs_ext->get_size() != ndof) {
    rhs_ext->free_data();
    rhs_ext->init(ndof);
  }
  else rhs_ext->set_zero();
  
  int k, m, marker;
  std::vector<AsmList> al(wf->neq);
  AsmList* am, * an;
  bool bnd[4];
  std::vector<bool> nat(wf->neq), isempty(wf->neq);
  EdgePos ep[4];
  reset_warn_order();

  if (rhsonly == false) {
    trace("Creating matrix sparse structure...");
    mat_ext->free_data();
  }
  else trace("Reusing matrix sparse structure...");

  if (rhsonly == false) trace("Assembling stiffness matrix and rhs...");
  else trace("Assembling rhs only...");
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
    for (unsigned int i = 0; i < s->idx.size(); i++) s->fns[i] = pss[s->idx[i]];
    for (unsigned int i = 0; i < s->ext.size(); i++) s->ext[i]->set_quad_2d(&g_quad_2d_std);
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
        // FIXME: Do not retrieve assembly list again if the element has not changed.
        spaces[j]->get_element_assembly_list(e[i], &al[j]);

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
        scalar **local_stiffness_matrix = get_matrix_buffer(std::max(am->cnt, an->cnt));
        for (int i = 0; i < am->cnt; i++)
        {
          if (!tra && am->dof[i] < 0) continue;
          fv->set_active_shape(am->idx[i]);

          if (!sym) // unsymmetric block
          {
            for (int j = 0; j < an->cnt; j++) {
              fu->set_active_shape(an->idx[j]);
              // FIXME - the NULL on the following eval_forms is temporary, an array of solutions 
              // should be passed there.
              if (an->dof[j] < 0) {
                if (dir_ext != NULL) {
                  scalar val = eval_form(mfv, NULL, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
                  dir_ext->add(am->dof[i], val);
                } 
              }
              else if (rhsonly == false) {
                scalar val = eval_form(mfv, NULL, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
                local_stiffness_matrix[i][j] = val;
              }
            }
          }
          else // symmetric block
          {
            for (int j = 0; j < an->cnt; j++) {
              if (j < i && an->dof[j] >= 0) continue;
              fu->set_active_shape(an->idx[j]);
              // FIXME - the NULL on the following eval_forms is temporary, an array of solutions 
              // should be passed there.
              if (an->dof[j] < 0) {
                if (dir_ext != NULL) {
                  scalar val = eval_form(mfv, NULL, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
                  dir_ext->add(am->dof[i], val);
                }
              } 
              else if (rhsonly == false) {
                scalar val = eval_form(mfv, NULL, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
                local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
              }
            }
          }
        }

        // insert the local stiffness matrix into the global one
        if (rhsonly == false) {
          insert_block(mat_ext, local_stiffness_matrix, am->dof, an->dof, am->cnt, an->cnt);
        }

        // insert also the off-diagonal (anti-)symmetric block, if required
        if (tra)
        {
          if (mfv->sym < 0) chsgn(local_stiffness_matrix, am->cnt, an->cnt);
          transpose(local_stiffness_matrix, am->cnt, an->cnt);
          if (rhsonly == false) {
            insert_block(mat_ext, local_stiffness_matrix, an->dof, am->dof, an->cnt, am->cnt);
          }

          // we also need to take care of the RHS...
          for (int j = 0; j < am->cnt; j++) {
            if (am->dof[j] < 0) {
              for (int i = 0; i < an->cnt; i++) {
                if (an->dof[i] >= 0) {
                  if (dir_ext != NULL) dir_ext->add(an->dof[i], local_stiffness_matrix[i][j]);
                }
              }
            }
          }
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
          scalar val = eval_form(vfv, NULL, fv, &refmap[m]) * am->coef[i];
          rhs_ext->add(am->dof[i], val);
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

          scalar **local_stiffness_matrix = get_matrix_buffer(std::max(am->cnt, an->cnt));
          for (int i = 0; i < am->cnt; i++)
          {
            if (am->dof[i] < 0) continue;
            fv->set_active_shape(am->idx[i]);
            for (int j = 0; j < an->cnt; j++)
            {
              fu->set_active_shape(an->idx[j]);
              // FIXME - the NULL on the following eval_forms is temporary, an array of solutions 
              // should be passed there.
              if (an->dof[j] < 0) {
                if (dir_ext != NULL) {
                  scalar val = eval_form(mfs, NULL, fu, fv, &refmap[n], &refmap[m], &(ep[edge])) 
                               * an->coef[j] * am->coef[i];
                  dir_ext->add(am->dof[i], val);
                }
              }
              else if (rhsonly == false) {
                scalar val = eval_form(mfs, NULL, fu, fv, &refmap[n], &refmap[m], &(ep[edge])) 
                             * an->coef[j] * am->coef[i];
                local_stiffness_matrix[i][j] = val;
              } 
            }
          }
          if (rhsonly == false) {
            insert_block(mat_ext, local_stiffness_matrix, am->dof, an->dof, am->cnt, an->cnt);
          }
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
            scalar val = eval_form(vfs, NULL, fv, &refmap[m], &(ep[edge])) * am->coef[i];
            rhs_ext->add(am->dof[i], val);
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

  if (rhsonly == false) values_changed = true;
}

// assembling of the jacobian matrix and residual vector for nonlinear problems
void DiscreteProblem::assemble(Matrix* mat_ext, Vector* rhs_ext, bool rhsonly)
{
  // the vector dir is irrelevant for nonlinear problems
  this->assemble(mat_ext, NULL, rhs_ext, rhsonly);
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

bool DiscreteProblem::solve_matrix_problem(Matrix* mat, Vector* vec) 
{
  // check matrix size
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof = 0 in DiscreteProblem::solve().");
  if (ndof != mat->get_size())
    error("Matrix size does not match ndof in DiscreteProblem:solve().");
  if (ndof != vec->get_size())
    error("Vector size does not match ndof in DiscreteProblem:solve().");

  // FIXME: similar test should be done for the vector "vec" also, but we need
  // to access the information about its length.
  // ...

  // solve the matrix problem (and report time)
  TimePeriod cpu_time;
  bool flag = this->solver->solve(mat, vec);
  report_time("Matrix problem solved in %g s", cpu_time.tick().last());

  return flag;
}

bool DiscreteProblem::solve(Matrix* mat, Vector* rhs, Vector* vec)
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
  Vector* delta = new AVector(ndof);
  memcpy(delta->get_c_array(), rhs, sizeof(scalar) * ndof);
  bool flag = this->solve_matrix_problem(mat, delta);
  if (flag == false) return false;

  // add the result which is in "delta" to the previous 
  // solution vector which is in "vec"
  for (int i = 0; i < ndof; i++) vec->add(i, delta->get_c_array()[i]);
  delete delta;

  return true;
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

/* OLD CODE
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
      warn("index of component: %d\n", i);
      error("Wrong projection norm in DiscreteProblem::project_global().");
    }
  }

  // FIXME: enable other types of matrices and vectors.
  CooMatrix mat(this->get_num_dofs());
  CommonSolverSciPyUmfpack solver;
  Vector* dir = new AVector(this->get_num_dofs());
  Vector* rhs = new AVector(this->get_num_dofs());

  //assembling the projection matrix, dir vector and rhs
  DiscreteProblem::assemble(&mat, dir, rhs, false);
  // since this is a linear problem, subtract the dir vector from the right-hand side:
  for (int i=0; i < this->get_num_dofs(); i++) rhs->add(i, -dir->get(i));

  solver.solve(&mat, rhs);
  for (int i=0; i < this->wf->neq; i++) {
    target[i]->set_fe_solution(this->spaces[i], rhs);
  }

  // restoring original weak form
  wf = wf_orig;
  wf_seq = -1;
}
*/

// global orthogonal projection
void project_global(Tuple<Space *> spaces, Tuple<MeshFunction*> source, 
                    Tuple<Solution*> target, WeakForm *wf, bool is_complex)
{
  int n = source.size();

  // sanity checks
  if (n <= 0 || n > 10) error("Wrong number of projected functions in project_global().");
  for (int i = 0; i < n; i++) if(spaces[i] == NULL) error("this->spaces[%d] == NULL in project_global().", i);
  if (spaces.size() != n) error("Number of spaces must matchnumber of projected functions in project_global().");
  if (target.size() != n) error("Mismatched numbers of projected functions and solutions in project_global().");

  // this is needed since spaces may have their DOFs enumerated only locally
  // when they come here.
  int ndof = assign_dofs(spaces);

  // FIXME: enable other types of matrices and vectors.
  CooMatrix mat(ndof, is_complex);
  CommonSolverSciPyUmfpack solver;
  Vector* dir = new AVector(ndof, is_complex);
  Vector* rhs = new AVector(ndof, is_complex);

  //assembling the projection matrix, dir vector and rhs  
  DiscreteProblem dp(wf, spaces);
  bool rhsonly = false;
  dp.assemble(&mat, dir, rhs, rhsonly);
  // since this is a linear problem, subtract the dir vector from the right-hand side:
  for (int i=0; i < ndof; i++) rhs->add(i, -dir->get(i));

  solver.solve(&mat, rhs);
  for (int i=0; i < wf->neq; i++) {
    target[i]->set_fe_solution(spaces[i], rhs);
  }
}

// global orthogonal projection
void project_global(Tuple<Space *> spaces, Tuple<MeshFunction*> source, 
                    Tuple<Solution*> target, Tuple<int>proj_norms, bool is_complex)
{
  int n = spaces.size();  

  // define temporary projection weak form
  WeakForm wf(n);
  int found[100];
  for (int i = 0; i < 100; i++) found[i] = 0;
  for (int i = 0; i < n; i++) {
    int norm;
    if (proj_norms == Tuple<int>()) norm = 1;
    else norm = proj_norms[i];
    if (norm == 0) {
      found[i] = 1;
      wf.add_matrix_form(i, i, L2projection_biform<double, scalar>, L2projection_biform<Ord, Ord>);
      wf.add_vector_form(i, L2projection_liform<double, scalar>, L2projection_liform<Ord, Ord>,
                     H2D_ANY, source[i]);
    }
    if (norm == 1) {
      found[i] = 1;
      wf.add_matrix_form(i, i, H1projection_biform<double, scalar>, H1projection_biform<Ord, Ord>);
      wf.add_vector_form(i, H1projection_liform<double, scalar>, H1projection_liform<Ord, Ord>,
                     H2D_ANY, source[i]);
    }
    if (norm == 2) {
      found[i] = 1;
      wf.add_matrix_form(i, i, Hcurlprojection_biform<double, scalar>, Hcurlprojection_biform<Ord, Ord>);
      wf.add_vector_form(i, Hcurlprojection_liform<double, scalar>, Hcurlprojection_liform<Ord, Ord>,
                     H2D_ANY, source[i]);
    }
  }
  for (int i=0; i < n; i++) {
    if (found[i] == 0) {
      warn("index of component: %d\n", i);
      error("Wrong projection norm in project_global().");
    }
  }

  project_global(spaces, source, target, &wf, is_complex);
}

void project_global(Tuple<Space *> spaces, Tuple<MeshFunction*> source, Tuple<Solution*> target,
		    matrix_forms_tuple_t proj_biforms, vector_forms_tuple_t proj_liforms, bool is_complex)
{
  int n = spaces.size();
  matrix_forms_tuple_t::size_type n_biforms = proj_biforms.size();
  if (n_biforms != proj_liforms.size())
    error("Mismatched numbers of projection forms in project_global().");
  if (n_biforms > 0) {
    if (n != n_biforms)
      error("Mismatched numbers of projected functions and projection forms in project_global().");
  }
  else
    warn("project_global() expected %d user defined biform(s) & liform(s); ordinary H1 projection will be performed", n);

  // this is needed since spaces may have their DOFs enumerated only locally
  // when they come here.
  int ndof = assign_dofs(spaces);

  // define projection weak form
  WeakForm wf(n);
  for (int i = 0; i < n; i++) {
    if (n_biforms == 0) {
      wf.add_matrix_form(i, i, H1projection_biform<double, scalar>, H1projection_biform<Ord, Ord>);
      wf.add_vector_form(i, H1projection_liform<double, scalar>, H1projection_liform<Ord, Ord>,
                     H2D_ANY, source[i]);
    }
    else {
      wf.add_matrix_form(i, i, proj_biforms[i].first, proj_biforms[i].second);
      wf.add_vector_form(i, proj_liforms[i].first, proj_liforms[i].second,
                     H2D_ANY, source[i]);
    }
  }

  project_global(spaces, source, target, &wf, is_complex);
}

/// Global orthogonal projection of one scalar ExactFunction.
void project_global(Space *space, ExactFunction source, Solution* target, int proj_norm, bool is_complex)
{
  if (proj_norm != 0 && proj_norm != 1) error("Wrong norm used in orthogonal projection (scalar case).");
  Mesh *mesh = space->get_mesh();
  if (mesh == NULL) error("Mesh is NULL in project_global().");
  Solution sln;
  sln.set_exact(mesh, source);
  project_global(space, &sln, target, proj_norm, is_complex);
};

/// Global orthogonal projection of one scalar ExactFunction -- user specified projection bi/linear forms.
void project_global(Space *space, ExactFunction source, Solution* target,
                    std::pair<WeakForm::matrix_form_val_t, WeakForm::matrix_form_ord_t> proj_biform,
                    std::pair<WeakForm::vector_form_val_t, WeakForm::vector_form_ord_t> proj_liform,
                    bool is_complex)
{
  // todo: check that supplied forms take scalar valued functions
  Mesh *mesh = space->get_mesh();
  if (mesh == NULL) error("Mesh is NULL in project_global().");
  Solution sln;
  sln.set_exact(mesh, source);
  project_global(space, &sln, target, matrix_forms_tuple_t(proj_biform), vector_forms_tuple_t(proj_liform), is_complex);
};

/// Global orthogonal projection of one vector-valued ExactFunction.
void project_global(Space *space, ExactFunction2 source, Solution* target, bool is_complex)
{
  int proj_norm = 2; // Hcurl
  Mesh *mesh = space->get_mesh();
  if (mesh == NULL) error("Mesh is NULL in project_global().");
  Solution sln;
  sln.set_exact(mesh, source);
  project_global(space, &sln, target, proj_norm, is_complex);
};

/// Projection-based interpolation of an exact function. This is faster than the
/// global projection since no global matrix problem is solved.
void project_local(Space *space, ExactFunction exactfn, Mesh* mesh,
                   Solution* result, int proj_norm, bool is_complex)
{
  /// TODO
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

/* TEMPORARILY DISABLED
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
    this->assemble(mat, dir, rhs, false);
    
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

*/
