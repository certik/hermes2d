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

void WeakForm::add_matrix_form(int i, int j, jacform_val_t fn, jacform_ord_t ord, SymFlag sym, int area, Tuple<MeshFunction*>ext)
{
  if (i < 0 || i >= neq || j < 0 || j >= neq)
    error("Invalid equation number.");
  if (sym < -1 || sym > 1)
    error("\"sym\" must be -1, 0 or 1.");
  if (sym < 0 && i == j)
    error("Only off-diagonal forms can be antisymmetric.");
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");
  if (jfvol.size() > 100)
    warn("Large number of forms (> 100). Is this the intent?");

  JacFormVol form = { i, j, sym, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  jfvol.push_back(form);
  seq++;
}

// single equation case
void WeakForm::add_matrix_form(jacform_val_t fn, jacform_ord_t ord, SymFlag sym, int area, Tuple<MeshFunction*>ext)
{
  int i = 0, j = 0;

  // FIXME: the code below should be replaced with a call to the full function.
  if (sym < -1 || sym > 1)
    error("\"sym\" must be -1, 0 or 1.");
  if (sym < 0 && i == j)
    error("Only off-diagonal forms can be antisymmetric.");
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");
  if (jfvol.size() > 100)
    warn("Large number of forms (> 100). Is this the intent?");

  JacFormVol form = { i, j, sym, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  jfvol.push_back(form);
  seq++;
}

void WeakForm::add_matrix_form_surf(int i, int j, jacform_val_t fn, jacform_ord_t ord, int area, Tuple<MeshFunction*>ext)
{
  if (i < 0 || i >= neq || j < 0 || j >= neq)
    error("Invalid equation number.");
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");

  JacFormSurf form = { i, j, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  jfsurf.push_back(form);
  seq++;
}

// single equation case
void WeakForm::add_matrix_form_surf(jacform_val_t fn, jacform_ord_t ord, int area, Tuple<MeshFunction*>ext)
{
  int i = 0, j = 0;

  // FIXME: the code below should be replaced with a call to the full function. 
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");

  JacFormSurf form = { i, j, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  jfsurf.push_back(form);
  seq++;
}

void WeakForm::add_vector_form(int i, resform_val_t fn, resform_ord_t ord, int area, Tuple<MeshFunction*>ext)
{
  if (i < 0 || i >= neq)
    error("Invalid equation number.");
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");

  ResFormVol form = { i, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  rfvol.push_back(form);
  seq++;
}

// single equation case
void WeakForm::add_vector_form(resform_val_t fn, resform_ord_t ord, int area, Tuple<MeshFunction*>ext)
{
  int i = 0;

  // FIXME: the code below should be replaced with a call to the full function. 
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");

  ResFormVol form = { i, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  rfvol.push_back(form);
  seq++;
}

void WeakForm::add_vector_form_surf(int i, resform_val_t fn, resform_ord_t ord, int area, Tuple<MeshFunction*>ext)
{
  if (i < 0 || i >= neq)
    error("Invalid equation number.");
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");

  ResFormSurf form = { i, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  rfsurf.push_back(form);
  seq++;
}

// single equation case
void WeakForm::add_vector_form_surf(resform_val_t fn, resform_ord_t ord, int area, Tuple<MeshFunction*>ext)
{
  int i = 0;

  // FIXME: the code below should be replaced with a call to the full function. 
  if (area != H2D_ANY && area < 0 && -area > areas.size())
    error("Invalid area number.");

  ResFormSurf form = { i, area, fn, ord };
  if (ext.size() != 0) {
    int nx = ext.size(); 
    for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
  }
  rfsurf.push_back(form);
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
void WeakForm::get_stages(Space** spaces, std::vector<WeakForm::Stage>& stages, bool rhsonly)
{
  unsigned i;
  stages.clear();

  // process volume jacforms
  for (i = 0; i < jfvol.size(); i++)
  {
    int ii = jfvol[i].i, jj = jfvol[i].j;
    Mesh* m1 = spaces[ii]->get_mesh();
    Mesh* m2 = spaces[jj]->get_mesh();
    Stage* s = find_stage(stages, ii, jj, m1, m2, jfvol[i].ext);
    s->jfvol.push_back(&jfvol[i]);
  }

  // process surface jacforms
  for (i = 0; i < jfsurf.size(); i++)
  {
    int ii = jfsurf[i].i, jj = jfsurf[i].j;
    Mesh* m1 = spaces[ii]->get_mesh();
    Mesh* m2 = spaces[jj]->get_mesh();
    Stage* s = find_stage(stages, ii, jj, m1, m2, jfsurf[i].ext);
    s->jfsurf.push_back(&jfsurf[i]);
  }

  // process volume res forms
  for (unsigned i = 0; i < rfvol.size(); i++) {
    int ii = rfvol[i].i;
    Mesh *m = spaces[ii]->get_mesh();
    Stage *s = find_stage(stages, ii, ii, m, m, rfvol[i].ext);
    s->rfvol.push_back(&rfvol[i]);
  }

  // process surface res forms
  for (unsigned i = 0; i < rfsurf.size(); i++) {
    int ii = rfsurf[i].i;
    Mesh *m = spaces[ii]->get_mesh();
    Stage *s = find_stage(stages, ii, ii, m, m, rfsurf[i].ext);
    s->rfsurf.push_back(&rfsurf[i]);
  }

  // helper macro for iterating in a set
  #define set_for_each(myset, type) \
    for (std::set<type>::iterator it = (myset).begin(); it != (myset).end(); it++)

  // initialize the arrays meshes and fns needed by Traverse for each stage
  for (i = 0; i < stages.size(); i++)
  {
    Stage* s = &stages[i];
    set_for_each(s->idx_set, int)
    {
      s->idx.push_back(*it);
      s->meshes.push_back(spaces[*it]->get_mesh());
      s->fns.push_back(NULL);
    }
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


/// Finds an assembling stage with the same set of meshes as [m1, m2, ext]. If no such
/// stage can be found, a new one is created and returned.
///
WeakForm::Stage* WeakForm::find_stage(std::vector<WeakForm::Stage>& stages, int ii, int jj,
                                      Mesh* m1, Mesh* m2, std::vector<MeshFunction*>& ext)
{
  // first create a list of meshes the form uses
  std::set<unsigned> seq;
  seq.insert(m1->get_seq());
  seq.insert(m2->get_seq());
  for (unsigned i = 0; i < ext.size(); i++)
    seq.insert(ext[i]->get_mesh()->get_seq());

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

  for (unsigned i = 0; i < jfvol.size(); i++) {
    blocks[jfvol[i].i][jfvol[i].j] = true;
    if (jfvol[i].sym)
      blocks[jfvol[i].j][jfvol[i].i] = true;
  }

  for (unsigned i = 0; i < jfsurf.size(); i++)
    blocks[jfsurf[i].i][jfsurf[i].j] = true;

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
