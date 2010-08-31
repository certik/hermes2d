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
#include "weakform.h"
#include "matrix_old.h"
#include "solution.h"
#include "forms.h"

//// interface /////////////////////////////////////////////////////////////////////////////////////

WeakForm::WeakForm(int neq, bool mat_free)
{
  this->neq = neq;
  seq = 0;
  this->is_matfree = mat_free;
}

void WeakForm::add_matrix_form(int i, int j, matrix_form_val_t fn, 
                               matrix_form_ord_t ord, SymFlag sym, int area, Tuple<MeshFunction*>ext)
{
  if (i < 0 || i >= neq || j < 0 || j >= neq)
    error("Invalid equation number.");
  if (sym < -1 || sym > 1)
    error("\"sym\" must be -1, 0 or 1.");
  if (sym < 0 && i == j)
    error("Only off-diagonal forms can be antisymmetric.");
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");
  if (mfvol.size() > 100)
    warn("Large number of forms (> 100). Is this the intent?");

  MatrixFormVol form = { i, j, sym, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  mfvol.push_back(form);
  seq++;
}

// single equation case
void WeakForm::add_matrix_form(matrix_form_val_t fn, matrix_form_ord_t ord, SymFlag sym, int area, Tuple<MeshFunction*>ext)
{
  int i = 0, j = 0;

  // FIXME: the code below should be replaced with a call to the full function.
  if (sym < -1 || sym > 1)
    error("\"sym\" must be -1, 0 or 1.");
  if (sym < 0 && i == j)
    error("Only off-diagonal forms can be antisymmetric.");
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");
  if (mfvol.size() > 100)
    warn("Large number of forms (> 100). Is this the intent?");

  MatrixFormVol form = { i, j, sym, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  mfvol.push_back(form);
  seq++;
}

void WeakForm::add_matrix_form_surf(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord, int area, Tuple<MeshFunction*>ext)
{
  if (i < 0 || i >= neq || j < 0 || j >= neq)
    error("Invalid equation number.");
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");

  MatrixFormSurf form = { i, j, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  mfsurf.push_back(form);
  seq++;
}

// single equation case
void WeakForm::add_matrix_form_surf(matrix_form_val_t fn, matrix_form_ord_t ord, int area, Tuple<MeshFunction*>ext)
{
  int i = 0, j = 0;

  // FIXME: the code below should be replaced with a call to the full function. 
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");

  MatrixFormSurf form = { i, j, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  mfsurf.push_back(form);
  seq++;
}

void WeakForm::add_vector_form(int i, vector_form_val_t fn, vector_form_ord_t ord, int area, Tuple<MeshFunction*>ext)
{
  if (i < 0 || i >= neq)
    error("Invalid equation number.");
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");

  VectorFormVol form = { i, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  vfvol.push_back(form);
  seq++;
}

// single equation case
void WeakForm::add_vector_form(vector_form_val_t fn, vector_form_ord_t ord, int area, Tuple<MeshFunction*>ext)
{
  int i = 0;

  // FIXME: the code below should be replaced with a call to the full function. 
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");

  VectorFormVol form = { i, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  vfvol.push_back(form);
  seq++;
}

void WeakForm::add_vector_form_surf(int i, vector_form_val_t fn, vector_form_ord_t ord, int area, Tuple<MeshFunction*>ext)
{
  if (i < 0 || i >= neq)
    error("Invalid equation number.");
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");

  VectorFormSurf form = { i, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  vfsurf.push_back(form);
  seq++;
}

// single equation case
void WeakForm::add_vector_form_surf(vector_form_val_t fn, vector_form_ord_t ord, int area, Tuple<MeshFunction*>ext)
{
  int i = 0;

  // FIXME: the code below should be replaced with a call to the full function. 
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");

  VectorFormSurf form = { i, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  vfsurf.push_back(form);
  seq++;
}

void WeakForm::set_ext_fns(void* fn, Tuple<MeshFunction*>ext)
{
  error("Not implemented yet.");
}


//// stages ////////////////////////////////////////////////////////////////////////////////////////

/// Constructs a list of assembling stages. Each stage contains a list of forms
/// that share the same meshes. Each stage is then assembled separately. This
/// improves the performance of multi-mesh assembling.
///
void WeakForm::get_stages(Tuple<Space *> spaces, Tuple<Solution *> u_ext, std::vector<WeakForm::Stage>& stages, bool rhsonly)
{
  unsigned i;
  stages.clear();

  // process volume matrix_forms
  for (i = 0; i < mfvol.size(); i++)
  {
    int ii = mfvol[i].i, jj = mfvol[i].j;
    Mesh* m1 = spaces[ii]->get_mesh();
    Mesh* m2 = spaces[jj]->get_mesh();
    Stage* s = find_stage(stages, ii, jj, m1, m2, mfvol[i].ext, u_ext);
    s->mfvol.push_back(&mfvol[i]);
  }

  // process surface matrix_forms
  for (i = 0; i < mfsurf.size(); i++)
  {
    int ii = mfsurf[i].i, jj = mfsurf[i].j;
    Mesh* m1 = spaces[ii]->get_mesh();
    Mesh* m2 = spaces[jj]->get_mesh();
    Stage* s = find_stage(stages, ii, jj, m1, m2, mfsurf[i].ext, u_ext);
    s->mfsurf.push_back(&mfsurf[i]);
  }

  // process volume res forms
  for (unsigned i = 0; i < vfvol.size(); i++) {
    int ii = vfvol[i].i;
    Mesh *m = spaces[ii]->get_mesh();
    Stage *s = find_stage(stages, ii, ii, m, m, vfvol[i].ext, u_ext);
    s->vfvol.push_back(&vfvol[i]);
  }

  // process surface res forms
  for (unsigned i = 0; i < vfsurf.size(); i++) {
    int ii = vfsurf[i].i;
    Mesh *m = spaces[ii]->get_mesh();
    Stage *s = find_stage(stages, ii, ii, m, m, vfsurf[i].ext, u_ext);
    s->vfsurf.push_back(&vfsurf[i]);
  }

  // helper macro for iterating in a set
  #define set_for_each(myset, type) \
    for (std::set<type>::iterator it = (myset).begin(); it != (myset).end(); it++)

  // initialize the arrays meshes and fns needed by Traverse for each stage
  for (i = 0; i < stages.size(); i++)
  {
    Stage* s = &stages[i];
    
    // First, initialize arrays for the test functions. A pointer to the PrecalcShapeset 
    // corresponding to each space will be assigned to s->fns later during assembling.
    set_for_each(s->idx_set, int)
    {
      s->idx.push_back(*it);
      s->meshes.push_back(spaces[*it]->get_mesh());
      s->fns.push_back(NULL);
    }
    
    // Next, append to the existing arrays the external functions (including the solutions
    // from previous Newton iteration) and their meshes. Also fill in a special array with
    // these external functions only.
    set_for_each(s->ext_set, MeshFunction*)
    {
      s->ext.push_back(*it);
      s->meshes.push_back((*it)->get_mesh());
      s->fns.push_back(*it);
    }
    
    s->idx_set.clear();
    s->seq_set.clear();
    s->ext_set.clear();
  }
}


/// Finds an assembling stage with the same set of meshes as [m1, m2, ext, u_ext]. If no such
/// stage can be found, a new one is created and returned.
///
WeakForm::Stage* WeakForm::find_stage(std::vector<WeakForm::Stage>& stages, int ii, int jj,
                                      Mesh* m1, Mesh* m2, 
                                      std::vector<MeshFunction*>& ext, std::vector<Solution*>& u_ext)
{
  // first create a list of meshes the form uses
  std::set<unsigned> seq;
  seq.insert(m1->get_seq());
  seq.insert(m2->get_seq());
  Mesh *mmm;
  for (unsigned i = 0; i < ext.size(); i++) {
    mmm = ext[i]->get_mesh();
    if (mmm == NULL) error("NULL Mesh pointer detected in ExtData during assembling.\n  Have you initialized all external functions?");
    seq.insert(mmm->get_seq());
  }
  for (unsigned i = 0; i < u_ext.size(); i++) {
    if (u_ext[i] != NULL) {
      mmm = u_ext[i]->get_mesh();
      if (mmm == NULL) error("NULL Mesh pointer detected in u_ext during assembling.");
      seq.insert(mmm->get_seq());
    }
  }
  
  // find a suitable existing stage for the form
  Stage* s = NULL;
  for (unsigned i = 0; i < stages.size(); i++)
    if (seq.size() == stages[i].seq_set.size() &&
        equal(seq.begin(), seq.end(), stages[i].seq_set.begin()))
      { s = &stages[i]; break; }

  // create a new stage if not found
  if (s == NULL)
  {
    Stage newstage;
    stages.push_back(newstage);
    s = &stages.back();
    s->seq_set = seq;
  }

  // update and return the stage
  for (unsigned int i = 0; i < ext.size(); i++)
    s->ext_set.insert(ext[i]);
  for (unsigned int i = 0; i < u_ext.size(); i++)
    if (u_ext[i] != NULL)
      s->ext_set.insert(u_ext[i]);
  
  s->idx_set.insert(ii);
  s->idx_set.insert(jj);
  return s;
}


/// Returns a (neq x neq) array containing true in each element, if the corresponding
/// block of weak forms is used, and false otherwise.
///
bool** WeakForm::get_blocks()
{
  bool** blocks = new_matrix<bool>(neq, neq);
  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
      blocks[i][j] = false;

  for (unsigned i = 0; i < mfvol.size(); i++) {
    blocks[mfvol[i].i][mfvol[i].j] = true;
    if (mfvol[i].sym)
      blocks[mfvol[i].j][mfvol[i].i] = true;
  }

  for (unsigned i = 0; i < mfsurf.size(); i++)
    blocks[mfsurf[i].i][mfsurf[i].j] = true;

  return blocks;
}


//// areas /////////////////////////////////////////////////////////////////////////////////////////

bool WeakForm::is_in_area_2(int marker, int area) const
{
  if (-area > (int)areas.size()) error("Invalid area number.");
  const Area* a = &areas[-area-1];

  for (unsigned int i = 0; i < a->markers.size(); i++)
    if (a->markers[i] == marker)
      return true;

  return false;
}
