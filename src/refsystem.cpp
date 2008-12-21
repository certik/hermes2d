// This file is part of Hermes2D.
//
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2005-2008 Pavel Solin <solin@utep.edu>
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
#include "space.h"
#include "weakform.h"
#include "refsystem.h"


  
RefSystem::RefSystem(LinSystem* coarse, int order_increase)
         : LinSystem(coarse->wf, coarse->solver)
{
  this->coarse = coarse;
  order_inc = order_increase;
  ref_meshes = NULL;
  ref_spaces = NULL;
}

RefSystem::~RefSystem()
{
  free_ref_data();
}


void RefSystem::set_spaces(int n, ...)
{
  error("set_spaces must not be called in RefSystem.");
}

void RefSystem::set_pss(int n, ...)
{
  error("set_pss must not be called in RefSystem.");
}


void RefSystem::assemble(bool rhsonly)
{
  int i, j;
  
  // get rid of any previous data
  free_ref_data();
  
  ref_meshes = new Mesh*[wf->neq];
  ref_spaces = new Space*[wf->neq];
  
  // copy meshes from the coarse problem, refine them
  for (i = 0; i < wf->neq; i++)
  {
    Mesh* cmesh = coarse->spaces[i]->get_mesh();
    
    // check if we already have the same mesh
    for (j = 0; j < i; j++)
      if (cmesh->get_seq() == coarse->spaces[j]->get_mesh()->get_seq())
        break;
      
    if (j < i) // yes
    {
      ref_meshes[i] = ref_meshes[j];
    }
    else // no, copy and refine the coarse one
    {
      Mesh* rmesh = new Mesh;
      rmesh->copy(cmesh);
      rmesh->refine_all_elements();
      ref_meshes[i] = rmesh;
    }
  }
  
  // duplicate spaces from the coarse problem, assign reference orders and dofs
  int dofs = 0;
  for (i = 0; i < wf->neq; i++)
  {
    ref_spaces[i] = coarse->spaces[i]->dup(ref_meshes[i]);
    ref_spaces[i]->copy_orders(coarse->spaces[i], order_inc);
    dofs += ref_spaces[i]->assign_dofs(dofs);
  }
  
  memcpy(spaces, ref_spaces, sizeof(Space*) * wf->neq);
  memcpy(pss, coarse->pss, sizeof(PrecalcShapeset*) * wf->neq);

  LinSystem::assemble(rhsonly);
}


void RefSystem::free_ref_data()
{
  int i, j;
  
  // free reference meshes
  if (ref_meshes != NULL)
  {
    for (i = 0; i < wf->neq; i++)
    {
      for (j = 0; j < i; j++)
        if (ref_meshes[j] == ref_meshes[i])
          break;
      
      if (i == j) delete ref_meshes[i];
    }
    
    delete [] ref_meshes;
    ref_meshes = NULL;
  }
  
  // free reference spaces
  if (ref_spaces != NULL)
  {
    for (i = 0; i < wf->neq; i++)
      delete ref_spaces[i];
    
    delete [] ref_spaces;
    ref_spaces = NULL;
  }
}
