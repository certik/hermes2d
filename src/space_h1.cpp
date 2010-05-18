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
#include "space_h1.h"
#include "matrix_old.h"
#include "quad_all.h"


double** H1Space::h1_proj_mat = NULL;
double*  H1Space::h1_chol_p   = NULL;
int      H1Space::h1_proj_ref = 0;


H1Space::H1Space(Mesh* mesh, Shapeset* shapeset)
       : Space(mesh, shapeset)
{
  if (!h1_proj_ref++)
  {
    // fixme: separate projection matrices for different shapesets
    precalculate_projection_matrix(2, h1_proj_mat, h1_chol_p);
  }
  proj_mat = h1_proj_mat;
  chol_p   = h1_chol_p;
}


H1Space::~H1Space()
{
  if (!--h1_proj_ref)
  {
    delete [] h1_proj_mat;
    delete [] h1_chol_p;
  }
}


Space* H1Space::dup(Mesh* mesh) const
{
  H1Space* space = new H1Space(mesh, shapeset);
  space->copy_callbacks(this);
  return space;
}


//// dof assignment ////////////////////////////////////////////////////////////////////////////////

void H1Space::assign_vertex_dofs()
{
  // Before assigning vertex DOFs, we must know which boundary vertex nodes are part of
  // a natural BC and which are part of an essential BC. The critical are those which
  // lie at an interface of both types of BCs and which must be treated as belonging
  // to the essential part. Unfortunately this has to be done on a per-space basis, as
  // the markers may have different meanings in different spaces. There is no way to
  // look at the adjacent edge nodes given a vertex node, thus we have to walk through
  // all elements in the mesh.

  // loop through all elements and assign vertex, edge and bubble dofs
  Element* e;
  for_all_active_elements(e, mesh)
  {
    int order = get_element_order(e->id);
    if (order > 0)
    {
      for (unsigned int i = 0; i < e->nvert; i++)
      {
        // vertex dofs
        Node* vn = e->vn[i];
        NodeData* nd = ndata + vn->id;
        if (!vn->is_constrained_vertex() && nd->dof == H2D_UNASSIGNED_DOF)
        {
          if (nd->n == BC_ESSENTIAL || is_fixed_vertex(vn->id))
          {
            nd->dof = H2D_CONSTRAINED_DOF;
          }
          else
          {
            nd->dof = next_dof;
            next_dof += stride;
          }
          nd->n = 1;
        }

        // edge dofs
        Node* en = e->en[i];
        nd = ndata + en->id;
        if (nd->dof == H2D_UNASSIGNED_DOF)
        {
          // if the edge node is not constrained, assign it dofs
          if (en->ref > 1 || en->bnd || mesh->peek_vertex_node(en->p1, en->p2) != NULL)
          {
            int ndofs = get_edge_order_internal(en) - 1;
            nd->n = ndofs;

            if (en->bnd && bc_type_callback(en->marker) == BC_ESSENTIAL)
            {
              nd->dof = H2D_CONSTRAINED_DOF;
            }
            else
            {
              nd->dof = next_dof;
              next_dof += ndofs * stride;
            }
          }
          else // constrained edge node
          {
            nd->n = -1;
          }
        }
      }
    }

    // bubble dofs
    shapeset->set_mode(e->get_mode());
    ElementData* ed = &edata[e->id];
    ed->bdof = next_dof;
    ed->n = order ? shapeset->get_num_bubbles(ed->order) : 0;
    next_dof += ed->n * stride;
  }
}


void H1Space::assign_edge_dofs()
{
  // TODO: remove this fn, we now assign all dofs at once
}

void H1Space::assign_bubble_dofs()
{
  // TODO: remove this fn, we now assign all dofs at once
}


//// assembly lists ////////////////////////////////////////////////////////////////////////////////

void H1Space::get_vertex_assembly_list(Element* e, int iv, AsmList* al)
{
  Node* vn = e->vn[iv];
  NodeData* nd = &ndata[vn->id];
  int index = shapeset->get_vertex_index(iv);
  if (get_element_order(e->id) == 0) return;

  if (!vn->is_constrained_vertex()) // unconstrained
  {
    al->add_triplet(index, nd->dof, (nd->dof >= 0) ? 1.0 : *(nd->vertex_bc_coef));
  }
  else // constrained
  {
    //debug_log("! B cause of the triplet\n");
    for (int j = 0; j < nd->ncomponents; j++)
      if (nd->baselist[j].coef != (scalar) 0)
      {
        al->add_triplet(index, nd->baselist[j].dof, nd->baselist[j].coef);
      }
  }
}


void H1Space::get_edge_assembly_list_internal(Element* e, int ie, AsmList* al)
{
  Node* en = e->en[ie];
  NodeData* nd = &ndata[en->id];
  if (get_element_order(e->id) == 0) return;

  if (nd->n >= 0) // unconstrained
  {
    if (nd->dof >= 0)
    {
      int ori = (e->vn[ie]->id < e->vn[e->next_vert(ie)]->id) ? 0 : 1;
      for (int j = 0, dof = nd->dof; j < nd->n; j++, dof += stride)
        al->add_triplet(shapeset->get_edge_index(ie, ori, j+2), dof, 1.0);
    }
    else
    {
      for (int j = 0; j < nd->n; j++)
      {
        al->add_triplet(shapeset->get_edge_index(ie, 0, j+2), -1, nd->edge_bc_proj[j+2]);
      }
    }
  }
  else // constrained
  {
    int part = nd->part;
    int ori = part < 0 ? 1 : 0;
    if (part < 0) part ^= ~0;

    nd = &ndata[nd->base->id];
    for (int j = 0, dof = nd->dof; j < nd->n; j++, dof += stride)
      al->add_triplet(shapeset->get_constrained_edge_index(ie, j+2, ori, part), dof, 1.0);
  }
}


//// BC stuff //////////////////////////////////////////////////////////////////////////////////////

scalar* H1Space::get_bc_projection(EdgePos* ep, int order)
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


//// hanging nodes /////////////////////////////////////////////////////////////////////////////////

inline void H1Space::output_component(BaseComponent*& current, BaseComponent*& last, BaseComponent* min,
                                      Node*& edge, BaseComponent*& edge_dofs)
{
  // if the dof is already in the list, just add half of the other coef
  if (last != NULL && last->dof == min->dof)
  {
    last->coef += min->coef * 0.5;
    return;
  }

  // leave space for edge node dofs if they belong in front of the current minimum dof
  if (edge != NULL && ndata[edge->id].dof <= min->dof)
  {
    edge_dofs = current;

    // (reserve space only if the edge dofs are not in the list yet)
    if (ndata[edge->id].dof != min->dof) {
      current += ndata[edge->id].n;
    }
    edge = NULL;
  }

  // output new dof
  current->dof = min->dof;
  current->coef = min->coef * 0.5;
  last = current++;
}


/// This is a documentation for sadfs
///
///
Space::BaseComponent* H1Space::merge_baselists(BaseComponent* l1, int n1, BaseComponent* l2, int n2,
                                               Node* edge, BaseComponent*& edge_dofs, int& ncomponents)
{
  // estimate the upper bound of the result size
  int max_result = n1 + n2;
  if (edge != NULL) max_result += ndata[edge->id].n;

  BaseComponent* result = (BaseComponent*) malloc(max_result * sizeof(BaseComponent));
  BaseComponent* current = result;
  BaseComponent* last = NULL;

  // main loop - always output the component with smaller dof so that we get a sorted array
  int i1 = 0, i2 = 0;
  while (i1 < n1 && i2 < n2)
  {
    if (l1[i1].dof < l2[i2].dof)
      output_component(current, last, l1 + i1++, edge, edge_dofs);
    else
      output_component(current, last, l2 + i2++, edge, edge_dofs);
  }

  // finish the longer baselist
  while (i1 < n1) output_component(current, last, l1 + i1++, edge, edge_dofs);
  while (i2 < n2) output_component(current, last, l2 + i2++, edge, edge_dofs);

  // don't forget to reserve space for edge dofs if we haven't done that already
  if (edge != NULL)
  {
    edge_dofs = current;
    current += ndata[edge->id].n;
  }

  // if we produced less components than we expected, reallocate the resulting array
  // ...this should be OK as we are always shrinking the array so no copying should occur
  ncomponents = current - result;
  if (ncomponents < max_result)
  {
    BaseComponent* reallocated_result = (BaseComponent*) realloc(result, ncomponents * sizeof(BaseComponent));
    if (edge_dofs != NULL)
    {
      edge_dofs = reallocated_result + (edge_dofs - result);
    }
    return reallocated_result;
  }
  else
    return result;
}


static Node* get_mid_edge_vertex_node(Element* e, int i, int j)
{
  if (e->is_triangle()) return e->sons[3]->vn[e->prev_vert(i)];
  else if (e->sons[2] == NULL) return i == 1 ? e->sons[0]->vn[2] : i == 3 ? e->sons[0]->vn[3] : NULL;
  else if (e->sons[0] == NULL) return i == 0 ? e->sons[2]->vn[1] : i == 2 ? e->sons[2]->vn[2] : NULL;
  else return e->sons[i]->vn[j];
}


void H1Space::update_constrained_nodes(Element* e, EdgeInfo* ei0, EdgeInfo* ei1, EdgeInfo* ei2, EdgeInfo* ei3)
{
  int j, k;
  EdgeInfo* ei[4] = { ei0, ei1, ei2, ei3 };
  NodeData* nd;

  if (get_element_order(e->id) == 0) return;

  // on non-refined elements all we have to do is update edge nodes lying on constrained edges
  if (e->active)
  {
    for (unsigned int i = 0; i < e->nvert; i++)
    {
      if (ei[i] != NULL)
      {
        nd = &ndata[e->en[i]->id];
        nd->base = ei[i]->node;
        nd->part = ei[i]->part;
        if (ei[i]->ori) nd->part ^= ~0;
      }
    }
  }
  // the element has sons - update mid-edge constrained vertex nodes
  else
  {
    // create new edge infos where we don't have them yet
    EdgeInfo ei_data[4];
    for (unsigned int i = 0; i < e->nvert; i++)
    {
      if (ei[i] == NULL)
      {
        j = e->next_vert(i);
        Node* mid_vn = get_mid_edge_vertex_node(e, i, j);
        if (mid_vn != NULL && mid_vn->is_constrained_vertex())
        {
          Node* mid_en = mesh->peek_edge_node(e->vn[i]->id, e->vn[j]->id);
          if (mid_en != NULL)
          {
            ei[i] = ei_data + i;
            ei[i]->node = mid_en;
            ei[i]->part = -1;
            ei[i]->lo = -1.0;
            ei[i]->hi =  1.0;
            ei[i]->ori = (e->vn[i]->id < e->vn[j]->id) ? 0 : 1;
          }
        }
      }
    }

    // create a baselist for each mid-edge vertex node
    for (unsigned int i = 0; i < e->nvert; i++)
    {
      if (ei[i] == NULL) continue;
      j = e->next_vert(i);

      Node* mid_vn = get_mid_edge_vertex_node(e, i, j);
      if (mid_vn == NULL) continue;

      Node* vn[2] = { e->vn[i], e->vn[j] }; // endpoint vertex nodes
      Node* en = ei[i]->node; // constraining edge node
      BaseComponent *bl[2], dummy_bl[2]; // base lists of v[0] and v[1]
      int nc[2] = { 0, 0 }; // number of components of bl[0] and bl[1]

      // get baselists of vn[0] and vn[1] - pretend we have them even if they are unconstrained
      for (k = 0; k < 2; k++)
      {
        nd = &ndata[vn[k]->id];
        if (vn[k]->is_constrained_vertex())
        {
          bl[k] = nd->baselist;
          nc[k] = nd->ncomponents;
        }
        else // make up an artificial baselist
        {
          dummy_bl[k].dof = nd->dof;
          dummy_bl[k].coef = (nd->dof >= 0) ? 1.0 : *nd->vertex_bc_coef;
          bl[k] = &dummy_bl[k];
          nc[k] = 1;
        }
      }

      // merge the baselists
      BaseComponent* edge_dofs;
      nd = &ndata[mid_vn->id];
      nd->baselist = merge_baselists(bl[0], nc[0], bl[1], nc[1], en, edge_dofs, nd->ncomponents);
      extra_data.push_back(nd->baselist);

      // set edge node coefs to function values of the edge functions
      double mid = (ei[i]->lo + ei[i]->hi) * 0.5;
      nd = &ndata[en->id];
      for (k = 0; k < nd->n; k++, edge_dofs++)
      {
        edge_dofs->dof = nd->dof + k*stride;
        edge_dofs->coef = shapeset->get_fn_value(shapeset->get_edge_index(0, ei[i]->ori, k+2), mid, -1.0, 0);
      }

      //dump_baselist(ndata[mid_vn->id]);
    }

    // create edge infos for half-edges
    EdgeInfo  half_ei_data[4][2];
    EdgeInfo* half_ei[4][2];
    for (unsigned int i = 0; i < e->nvert; i++)
    {
      if (ei[i] == NULL)
      {
        half_ei[i][0] = half_ei[i][1] = NULL;
      }
      else
      {
        EdgeInfo* h0 = half_ei[i][0] = half_ei_data[i];
        EdgeInfo* h1 = half_ei[i][1] = half_ei_data[i] + 1;

        h0->node = h1->node = ei[i]->node;
        h0->part = (ei[i]->part + 1) * 2;
        h1->part = h0->part + 1;
        h0->hi = h1->lo = (ei[i]->lo + ei[i]->hi) / 2;
        h0->lo = ei[i]->lo;
        h1->hi = ei[i]->hi;
        h1->ori = h0->ori = ei[i]->ori;
      }
    }

    // recur to sons
    if (e->is_triangle())
    {
      update_constrained_nodes(e->sons[0], half_ei[0][0], NULL, half_ei[2][1], NULL);
      update_constrained_nodes(e->sons[1], half_ei[0][1], half_ei[1][0], NULL, NULL);
      update_constrained_nodes(e->sons[2], NULL, half_ei[1][1], half_ei[2][0], NULL);
      update_constrained_nodes(e->sons[3], NULL, NULL, NULL, NULL);
    }
    else if (e->sons[2] == NULL) // 'horizontally' split quad
    {
      update_constrained_nodes(e->sons[0], ei[0], half_ei[1][0], NULL, half_ei[3][1]);
      update_constrained_nodes(e->sons[1], NULL, half_ei[1][1], ei[2], half_ei[3][0]);
    }
    else if (e->sons[0] == NULL) // 'vertically' split quad
    {
      update_constrained_nodes(e->sons[2], half_ei[0][0], NULL, half_ei[2][1], ei[3]);
      update_constrained_nodes(e->sons[3], half_ei[0][1], ei[1], half_ei[2][0], NULL);
    }
    else // fully split quad
    {
      update_constrained_nodes(e->sons[0], half_ei[0][0], NULL, NULL, half_ei[3][1]);
      update_constrained_nodes(e->sons[1], half_ei[0][1], half_ei[1][0], NULL, NULL);
      update_constrained_nodes(e->sons[2], NULL, half_ei[1][1], half_ei[2][0], NULL);
      update_constrained_nodes(e->sons[3], NULL, NULL, half_ei[2][1], half_ei[3][0]);
    }
  }
}


void H1Space::update_constraints()
{
  Element* e;
  for_all_base_elements(e, mesh)
    update_constrained_nodes(e, NULL, NULL, NULL, NULL);
}


/*void H1Space::dump_baselist(NodeData& nd)
{
  printf("  { ");
  for (int i = 0; i < nd.ncomponents; i++)
    printf("{ %d, %lg } ", nd.baselist[i].dof, nd.baselist[i].coef);
  printf(" }\n");
}*/


//// vertex fixing /////////////////////////////////////////////////////////////////////////////////

void H1Space::fix_vertex(int id, scalar value)
{
  FixedVertex fv = { id, value };
  fixed_vertices.push_back(fv);
}


bool H1Space::is_fixed_vertex(int id) const
{
  for (unsigned int i = 0; i < fixed_vertices.size(); i++)
    if (fixed_vertices[i].id == id)
      return true;

  return false;
}


void H1Space::post_assign()
{
  // process fixed vertices -- put their values into nd->vertex_bc_coef
  for (unsigned int i = 0; i < fixed_vertices.size(); i++)
  {
    scalar* fixv = new scalar[1];
    *fixv = fixed_vertices[i].value;
    NodeData* nd = &ndata[fixed_vertices[i].id];
    nd->vertex_bc_coef = fixv;
    extra_data.push_back(fixv);
  }
}

