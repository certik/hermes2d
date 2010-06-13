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
  if (this->base->have_spaces == false) 
    error("RefSystem: spaces in base system not up to date.");

  this->wf = this->base->wf;
  this->order_increase = order_increase;
  this->refinement = refinement;
  this->linear = base->is_linear();

  // Have to set it manually as Nonlinsystem constructor
  // always sets it false;
  if (this->linear) {
    this->want_dir_contrib = true;
  } else {
    this->want_dir_contrib = false;
  }
}

void RefSystem::prepare() 
{
  // free meshes and spaces
  this->free_meshes_and_spaces();

  // create new meshes and spaces
  this->meshes = new Mesh*[this->wf->neq];
  this->spaces = new Space*[this->wf->neq];
  this->sp_seq = new int[this->wf->neq];
  for (int i=0; i < this->wf->neq; i++) {
    this->spaces[i]->set_seq(this->base->spaces[i]->get_seq());
    this->sp_seq[i] = this->base->spaces[i]->get_seq();
  }

  int i, j;
  // copy meshes from the coarse problem and refine them
  for (i = 0; i < this->wf->neq; i++)
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
      this->meshes[i] = rmesh;
    }
  }

  // duplicate spaces from the coarse problem, assign reference orders and dofs
  int ndof = 0;
  for (i = 0; i < this->wf->neq; i++)
  {
    this->spaces[i] = this->base->spaces[i]->dup(this->meshes[i]);

    if (refinement == -1)
    {
      Element* re;
      for_all_active_elements(re, meshes[i])
      {
        Mesh* mesh = this->base->spaces[i]->get_mesh();
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
    else {
      spaces[i]->copy_orders(base->spaces[i], order_increase);
    }

    ndof += spaces[i]->assign_dofs(ndof);
  }

  //memcpy(spaces, base->spaces, sizeof(Space*) * this->wf->neq);
  //memcpy(pss, base->pss, sizeof(PrecalcShapeset*) * this->wf->neq);
  
  this->pss = new PrecalcShapeset*[this->wf->neq];
  for(int i=0; i < this->wf->neq; i++) this->pss[i] = this->base->pss[i];
  this->have_spaces = true;
}

RefSystem::~RefSystem()
{
  // Nothing to do here, all is freed in the destructor to LinSystem.
}

// just to prevent the user from calling this
void RefSystem::set_spaces(Tuple<Space*> spaces)
{
  error("Method set_spaces must not be called in RefSystem.");
}

// just to prevent the user from calling this
void RefSystem::set_pss(Tuple<PrecalcShapeset*> pss)
{
  error("Method set_pss must not be called in RefSystem.");
}

void RefSystem::set_order_increase(int order_increase)
{
  this->order_increase = order_increase;
}

void RefSystem::assemble(bool rhsonly)
{  
  // perform uniform mesh refinement
  printf("prepare refsystem 1\n");
  prepare();
  printf("prepare refsystem 2\n");

  // internal check
  if (this->have_spaces == false) error("RefSystem: missing space(s).");

  // call LinSystem's assemble() function.
  if (!linear) {
    NonlinSystem::assemble(rhsonly);
  } else {
    LinSystem::assemble(rhsonly);
  }
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


