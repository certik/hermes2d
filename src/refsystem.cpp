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
#include "space.h"
#include "weakform.h"
#include "refsystem.h"
#include "solution.h"

RefSystem::RefSystem(LinSystem* base, int order_increase,
    int refinement) : NonlinSystem(base->wf, base->solver)
{
  this->base = base;
  this->spaces = NULL;
  this->num_spaces = base->get_num_spaces();
  this->order_increase = order_increase;
  this->refinement = refinement;
  this->nonlinear = !base->is_linear();
}


RefSystem::~RefSystem()
{
  // Nothing to do here, all is freed in the destructor to LinSystem.
}

void RefSystem::set_spaces(int n, ...)
{
  error("Method set_spaces must not be called in RefSystem.");
}

void RefSystem::set_pss(int n, ...)
{
  error("Method set_pss must not be called in RefSystem.");
}

void RefSystem::set_order_increase(int order_increase)
{
  this->order_increase = order_increase;
}

void RefSystem::assemble(bool rhsonly)
{  
  // copy and refines meshes from coarse mesh linsystem
  prepare();

  // call LinSystem's assemble() function.
  if (nonlinear) {
    NonlinSystem::assemble(rhsonly);
  } else {
    LinSystem::assemble(rhsonly);
  }
}

void RefSystem::prepare()
{
  int i, j;

  // free meshes and spaces
  free_meshes_and_spaces();

  // create new meshes and spaces
  meshes = new Mesh*[this->num_spaces];
  spaces = new Space*[this->num_spaces];

  // copy meshes from the coarse problem and refine them
  for (i = 0; i < this->num_spaces; i++)
  {
    Mesh* mesh = base->spaces[i]->get_mesh();

    // check if we already have the same mesh
    for (j = 0; j < i; j++)
      if (mesh->get_seq() == base->spaces[j]->get_mesh()->get_seq())
        break;

    if (j < i) // yes
    {
      meshes[i] = meshes[j];
    }
    else // no, copy and refine the coarse one
    {
      Mesh* rmesh = new Mesh;
      rmesh->copy(mesh);
      if (refinement ==  1) rmesh->refine_all_elements();
      if (refinement == -1) rmesh->unrefine_all_elements();
      meshes[i] = rmesh;
    }
  }

  // duplicate spaces from the coarse problem, assign reference orders and dofs
  int dofs = 0;
  for (i = 0; i < this->num_spaces; i++)
  {
    spaces[i] = base->spaces[i]->dup(meshes[i]);
    if (refinement == -1)
    {
      Element* re;
      for_all_active_elements(re, meshes[i])
      {
        Mesh* mesh = base->spaces[i]->get_mesh();
        Element* e = mesh->get_element(re->id);
        int max_order_h = 0, max_order_v = 0;
        if (e->active) {
          int quad_order = base->spaces[i]->get_element_order(e->id);
          max_order_h = H2D_GET_H_ORDER(quad_order);
          max_order_v = H2D_GET_V_ORDER(quad_order);
        }
        else
        { //find maximum order of sons
          for (int son = 0; son < 4; son++)
          {
            if (e->sons[son] != NULL)
            {
              int quad_order = base->spaces[i]->get_element_order(e->sons[son]->id);
              max_order_h = std::max(max_order_h, H2D_GET_H_ORDER(quad_order));
              max_order_v = std::max(max_order_v, H2D_GET_V_ORDER(quad_order));
            }
          }
        }

        //increase order and set it to element
        max_order_h = std::max(1, max_order_h + order_increase);
        if (re->is_triangle())
          max_order_v = 0;
        else
          max_order_v = std::max(1, max_order_v + order_increase);
        spaces[i]->set_element_order(re->id, H2D_MAKE_QUAD_ORDER(max_order_h, max_order_v));
      }
    }
    else
      spaces[i]->copy_orders(base->spaces[i], order_increase);

    dofs += spaces[i]->assign_dofs(dofs);
  }

  this->ndofs = dofs;
  memcpy(spaces, spaces, sizeof(Space*) * this->num_spaces);
  memcpy(pss, base->pss, sizeof(PrecalcShapeset*) * this->num_spaces);
  have_spaces = true;
}

bool RefSystem::solve_exact(scalar (*exactfn)(double x, double y, scalar& dx , scalar& dy), Solution* sln)
{
  Space* space = spaces[0];

  // some sanity checkz
  //if (!space->is_up_to_date())
  //  error("'space' is not up to date.");

  //set mesh and function
  sln->set_exact(meshes[0], exactfn);

  return true;
}


