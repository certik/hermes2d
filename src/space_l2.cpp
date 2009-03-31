// This file is part of Hermes2D.
//
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2005-2008 Pavel Solin <solin@unr.edu>
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

// $Id: space_l2.cpp 1086 2008-10-21 09:05:44Z jakub $

#include "common.h"
#include "space_l2.h"
#include "matrix.h"
#include "quad_all.h"


double** L2Space::l2_proj_mat = NULL;
double*  L2Space::l2_chol_p   = NULL;
int      L2Space::l2_proj_ref = 0;


L2Space::L2Space(Mesh* mesh, Shapeset* shapeset)
       : Space(mesh, shapeset)
{
  if (!l2_proj_ref++)
  {
    precalculate_projection_matrix(2, l2_proj_mat, l2_chol_p);
  }
  proj_mat = l2_proj_mat;
  chol_p   = l2_chol_p;

  ldata = NULL;
  lsize = 0;
}


L2Space::~L2Space()
{
  if (!--l2_proj_ref)
  {
    delete [] l2_proj_mat;
    delete [] l2_chol_p;
  }
  ::free(ldata);
}


Space* L2Space::dup(Mesh* mesh) const
{
  L2Space* space = new L2Space(mesh, shapeset);
  space->copy_callbacks(this);
  return space;
}


int L2Space::get_edge_order(Element* e, int edge)
{
  int order = edata[e->id].order;
  if (e->is_triangle()) return order;
  return (edge & 1) ? get_v_order(order) : get_h_order(order);
}


//// dof assignment ////////////////////////////////////////////////////////////////////////////////

void L2Space::resize_tables()
{
  if (lsize < mesh->get_max_element_id())
  {
    if (!lsize) lsize = 1000;
    while (lsize < mesh->get_max_element_id()) lsize = lsize * 3 / 2;
    ldata = (L2Data*) realloc(ldata, sizeof(L2Data) * lsize);
  }
  Space::resize_tables();
}


void L2Space::assign_vertex_dofs()
{
  int i, j;
  Element* e;
  for_all_active_elements(e, mesh)
  {
    // determine boundary types of this element
    int bc_type[4];
    for (i = 0; i < e->nvert; i++)
      if (!e->en[i]->bnd)
        bc_type[i] = BC_NONE;
      else
        bc_type[i] = bc_type_callback(e->en[i]->marker);
 
    // assign vertex dofs where appropriate
    L2Data* ld = &ldata[e->id];
    for (i = 0; i < e->nvert; i++)
    {
      j = e->prev_vert(i);
      if (bc_type[i] == BC_ESSENTIAL || bc_type[j] == BC_ESSENTIAL)
      {
        ld->vdof[i] = -1;
      }
      else
      {
        ld->vdof[i] = next_dof;
        next_dof += stride;
      }
    }
  }  
}


void L2Space::assign_edge_dofs()
{
  Element* e;
  for_all_active_elements(e, mesh)
  {
    // assign edge dofs on non-Dirichlet edges
    L2Data* ld = &ldata[e->id];
    for (int i = 0; i < e->nvert; i++)
    {
      int eo = edata[e->id].order;
      if (e->is_quad() && (i == 0 || i == 2))
        eo = get_v_order(eo);
      else
        eo = get_h_order(eo);

      Node* en = e->en[i];
      if (en->bnd && bc_type_callback(en->marker) == BC_ESSENTIAL)
      {
        ld->edof[i] = -1;
      }
      else
      {
        ld->edof[i] = next_dof;
        next_dof += (eo - 1) * stride;
      }
    }
  }
}


void L2Space::assign_bubble_dofs()
{
  Element* e;
  for_all_active_elements(e, mesh)
  {
    shapeset->set_mode(e->get_mode());
    ElementData* ed = &edata[e->id];
    ed->bdof = next_dof;
    ed->n = shapeset->get_num_bubbles(ed->order);
    next_dof += ed->n * stride;    
  }
}


//// assembly lists ////////////////////////////////////////////////////////////////////////////////

void L2Space::get_vertex_assembly_list(Element* e, int iv, AsmList* al)
{
  L2Data* ld = &ldata[e->id];
  
  scalar coef = 1.0;
  if (ld->vdof[iv] < 0)
  {
    /*if (bc_type_callback(e->en[iv]->marker) == BC_ESSENTIAL)
      coef = ndata[e->en[iv]->id].edge_bc_proj[0];
    else
      coef = ndata[e->en[e->prev_vert(iv)]->id].edge_bc_proj[1];*/
    coef = *(ndata[e->vn[iv]->id].vertex_bc_coef);
  }
  
  al->add_triplet(shapeset->get_vertex_index(iv), ld->vdof[iv], coef);
}


void L2Space::get_edge_assembly_list_internal(Element* e, int ie, AsmList* al)
{
  int eo = edata[e->id].order;
  if (e->is_quad() && (ie == 0 || ie == 2))
    eo = get_v_order(eo);
  else
    eo = get_h_order(eo);
  
  L2Data* ld = &ldata[e->id];
  if (ld->edof[ie] >= 0)
  {
    for (int j = 0, dof = ld->edof[ie]; j < eo-1; j++, dof += stride)
      al->add_triplet(shapeset->get_edge_index(ie, 0, j+2), dof, 1.0);
  }
  else
  {
    NodeData* nd = &ndata[e->en[ie]->id];
    for (int j = 0; j < eo-1; j++)
      al->add_triplet(shapeset->get_edge_index(ie, 0, j+2), -1, nd->edge_bc_proj[j+2]);
  }
}


void L2Space::get_bubble_assembly_list(Element* e, AsmList* al)
{
  ElementData* ed = &edata[e->id];
  if (!ed->n) return;

  int* indices = shapeset->get_bubble_indices(ed->order);
  for (int i = 0, dof = ed->bdof; i < ed->n; i++, dof += stride)
    al->add_triplet(*indices++, dof, 1.0);
}


//// BC stuff //////////////////////////////////////////////////////////////////////////////////////

scalar* L2Space::get_bc_projection(EdgePos* ep, int order)
{
  assert(order >= 1);
  scalar* proj = new scalar[order + 1];
  
  // obtain linear part of the projection
  ep->t = ep->lo;
  proj[0] = bc_value_callback_by_edge(ep);
  ep->t = ep->hi;
  proj[1] = bc_value_callback_by_edge(ep);

  if (order-- > 1)
  {
    Quad1DStd quad1d;
    scalar* rhs = proj + 2;
    int mo = quad1d.get_max_order();
    double2* pt = quad1d.get_points(mo);

    // get boundary values at integration points, construct rhs
    for (int i = 0; i < order; i++)
    {
      rhs[i] = 0.0;
      int ii = shapeset->get_edge_index(0, 0, i+2);
      for (int j = 0; j < quad1d.get_num_points(mo); j++)
      {
        double t = (pt[j][0] + 1) * 0.5, s = 1.0 - t;
        scalar l = proj[0] * s + proj[1] * t;
        ep->t = ep->lo * s + ep->hi * t;
        rhs[i] += pt[j][1] * shapeset->get_fn_value(ii, pt[j][0], -1.0, 0)
                           * (bc_value_callback_by_edge(ep) - l);
      }
    }

    // solve the system using a precalculated Cholesky decomposed projection matrix
    cholsl(proj_mat, order, chol_p, rhs, rhs);
  }

  return proj;
}
