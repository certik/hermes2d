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
#include "mesh.h"


//// nodes, element ////////////////////////////////////////////////////////////////////////////////

void Node::ref_element(Element* e)
{
  if (type == TYPE_EDGE)
  {
    // store the element pointer in a free slot of 'elem'
    if (elem[0] == NULL) elem[0] = e;
    else if (elem[1] == NULL) elem[1] = e;
    else assert(0);
  }
  ref++;
}


void Node::unref_element(HashTable* ht, Element* e)
{
  if (type == TYPE_VERTEX)
  {
    if (!--ref) ht->remove_vertex_node(id);
  }
  else
  {
    // remove the element from the array 'elem'
    if (elem[0] == e) elem[0] = NULL;
    else if (elem[1] == e) elem[1] = NULL;

    if (!--ref) ht->remove_edge_node(id);
  }
}


void Element::ref_all_nodes()
{
  for (int i = 0; i < nvert; i++)
  {
    vn[i]->ref_element();
    en[i]->ref_element(this);
  }
}


void Element::unref_all_nodes(HashTable* ht)
{
  for (int i = 0; i < nvert; i++)
  {
    vn[i]->unref_element(ht);
    en[i]->unref_element(ht, this);
  }
}


Element* Element::get_neighbor(int ie) const
{
  Element** elem = en[ie]->elem;
  if (elem[0] == this) return elem[1];
  if (elem[1] == this) return elem[0];
  assert(0);
}


double Element::get_area() const
{
  double ax, ay, bx, by;
  ax = vn[1]->x - vn[0]->x;
  ay = vn[1]->y - vn[0]->y;
  bx = vn[2]->x - vn[0]->x;
  by = vn[2]->y - vn[0]->y;

  double area = 0.5*(ax*by - ay*bx);
  if (is_triangle()) return area;

  ax = bx; ay = by;
  bx = vn[3]->x - vn[0]->x;
  by = vn[3]->y - vn[0]->y;

  return area + 0.5*(ax*by - ay*bx);
}


double Element::get_diameter() const
{
  double max, l;
  if (is_triangle())
  {
    max = 0.0;
    for (int i = 0; i < 3; i++)
    {
      int j = next_vert(i);
      l = sqr(vn[i]->x - vn[j]->x) + sqr(vn[i]->y - vn[j]->y);
      if (l > max) max = l;
    }
  }
  else
  {
    max = sqr(vn[0]->x - vn[2]->x) + sqr(vn[0]->y - vn[2]->y);
    l   = sqr(vn[1]->x - vn[3]->x) + sqr(vn[1]->y - vn[3]->y);
    if (l > max) max = l;
  }
  return sqrt(max);
}


//// mesh //////////////////////////////////////////////////////////////////////////////////////////

unsigned g_mesh_seq = 0;


Mesh::Mesh() : HashTable()
{
  nbase = nactive = ntopvert = ninitial = 0;
  seq = g_mesh_seq++;
}


Element* Mesh::get_element(int id) const
{
  if (id < 0 || id >= elements.get_size())
    error("invalid element id number.");
  return &(elements[id]);
}


int Mesh::get_edge_sons(Element* e, int edge, int& son1, int& son2)
{
  assert(!e->active);

  if (!e->is_triangle())
  {
    if (e->sons[2] == NULL) // horz quad
    {
      if (edge == 0 || edge == 2) { son1 = edge >> 1;   return 1; }
              else if (edge == 1) { son1 = 0; son2 = 1; return 2; }
                             else { son1 = 1; son2 = 0; return 2; }
    }
    else if (e->sons[0] == NULL) // vert quad
    {
      if (edge == 1 || edge == 3) { son1 = (edge == 1) ? 3 : 2; return 1; }
              else if (edge == 0) { son1 = 2; son2 = 3; return 2; }
                             else { son1 = 3; son2 = 2; return 2; }
    }
  }

  // triangle or 4-son quad
  son1 = edge;
  son2 = e->next_vert(edge);
  return 2;
}


//// low-level refinement //////////////////////////////////////////////////////////////////////////

/*  node and son numbering on a triangle:

                    vn[2]                                       vn[2]

                      *                                           *
                     / \                                         / \
                    /   \                                       /   \
                   /     \                                     /     \
                  /       \                                   / son[2]\
                 /         \                                 /_________\
        en[2]   /           \   en[1]                 vn[0] *           * vn[1]
               *             *                       vn[1]  *-----------*  vn[0]
              /               \                     vn[2] *  \         /  * vn[2]
             /                 \                         / \  \ son[3]/  / \
            /                   \                       /   \  \     /  /   \
           /                     \                     /     \  \   /  /     \
          /                       \                   / son[0]\  \ /  /son[1] \
         /                         \                 /         \  *  /         \
        *-------------*-------------*               *-----------*   *-----------*
                                               vn[0]      vn[1] vn[2] vn[0]      vn[1]
    vn[0]           en[0]           vn[1]


   node and son numbering on a quad:          refinement '0':

    vn[3]           en[2]           vn[2]       vn[3]        vn[2] vn[3]        vn[2]

        *-------------*-------------*               *------------* *------------*
        |                           |               |            | |            |
        |                           |               |            | |            |
        |                           |               |   son[3]   | |   son[2]   |
        |                           |               |            | |            |
        |                           |               |       vn[1]| |vn[0]       |
        |                           |         vn[0] *------------* *------------* vn[1]
 en[3]  *                           *  en[1]  vn[3] *------------* *------------* vn[2]
        |                           |               |       vn[2]| |vn[3]       |
        |                           |               |            | |            |
        |                           |               |   son[0]   | |   son[1]   |
        |                           |               |            | |            |
        |                           |               |            | |            |
        |                           |               *------------* *------------*
        *-------------*-------------*
                                                vn[0]        vn[1] vn[0]        vn[1]
    vn[0]           en[0]           vn[1]


  refinement '1':                             refinement '2':

    vn[3]                           vn[2]       vn[3]        vn[2] vn[3]        vn[2]

        *---------------------------*               *------------* *------------*
        |                           |               |            | |            |
        |                           |               |            | |            |
        |          son[1]           |               |            | |            |
        |                           |               |            | |            |
        |                           |               |            | |            |
  vn[0] *---------------------------* vn[1]         |            | |            |
  vn[3] *---------------------------* vn[2]         |   son[2]   | |   son[3]   |
        |                           |               |            | |            |
        |                           |               |            | |            |
        |          son[0]           |               |            | |            |
        |                           |               |            | |            |
        |                           |               |            | |            |
        *---------------------------*               *------------* *------------*

    vn[0]                           vn[1]       vn[0]        vn[1] vn[0]        vn[1]

*/

Element* Mesh::create_triangle(int marker, Node* v0, Node* v1, Node* v2, CurvMap* cm)
{
  // create a new element
  Element* e = elements.add();
  e->active = 1;
  e->marker = marker;
  e->userdata = 0;
  e->nvert = 3;
  e->iro_cache = -1;
  e->cm = cm;

  // set vertex and edge node pointers
  e->vn[0] = v0;
  e->vn[1] = v1;
  e->vn[2] = v2;
  e->en[0] = get_edge_node(v0->id, v1->id);
  e->en[1] = get_edge_node(v1->id, v2->id);
  e->en[2] = get_edge_node(v2->id, v0->id);

  // register in the nodes
  e->ref_all_nodes();

  return e;
}


Element* Mesh::create_quad(int marker, Node* v0, Node* v1, Node* v2, Node* v3, CurvMap* cm)
{
  // create a new element
  Element* e = elements.add();
  e->active = 1;
  e->marker = marker;
  e->userdata = 0;
  e->nvert = 4;
  e->iro_cache = -1;
  e->cm = cm;

  // set vertex and edge node pointers
  e->vn[0] = v0;
  e->vn[1] = v1;
  e->vn[2] = v2;
  e->vn[3] = v3;
  e->en[0] = get_edge_node(v0->id, v1->id);
  e->en[1] = get_edge_node(v1->id, v2->id);
  e->en[2] = get_edge_node(v2->id, v3->id);
  e->en[3] = get_edge_node(v3->id, v0->id);

  // register in the nodes
  e->ref_all_nodes();

  return e;
}


static CurvMap* create_son_curv_map(Element* e, int son)
{
  // if the top three bits of part are nonzero, we would overflow
  // -- make the element non-curvilinear
  if (e->cm->part & 0xe000000000000000ULL) return NULL;

  // if the parent element is already almost straight-edged,
  // the son will be even more straight-edged
  if (e->iro_cache == 0) return NULL;

  CurvMap* cm = new CurvMap;
  if (e->cm->toplevel == false)
  {
    cm->parent = e->cm->parent;
    cm->part = (e->cm->part << 3) + son + 1;
  }
  else
  {
    cm->parent = e;
    cm->part = (son + 1);
  }
  cm->toplevel = false;
  cm->order = 4;

  return cm;
}


void Mesh::refine_triangle(Element* e)
{
  // remember the markers of the edge nodes
  int bnd[3] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd    };
  int mrk[3] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker };

  // obtain three mid-edge vertex nodes
  Node* x0 = get_vertex_node(e->vn[0]->id, e->vn[1]->id);
  Node* x1 = get_vertex_node(e->vn[1]->id, e->vn[2]->id);
  Node* x2 = get_vertex_node(e->vn[2]->id, e->vn[0]->id);

  CurvMap* cm[4];
  memset(cm, 0, sizeof(cm));

  // adjust mid-edge coordinates if this is a curved element
  if (e->is_curved())
  {
    double2 pt[3] = { { 0.0,-1.0 }, { 0.0, 0.0 }, { -1.0, 0.0 } };
    e->cm->get_mid_edge_points(e, pt, 3);
    x0->x = pt[0][0]; x0->y = pt[0][1];
    x1->x = pt[1][0]; x1->y = pt[1][1];
    x2->x = pt[2][0]; x2->y = pt[2][1];

    // create CurvMaps for sons (pointer to parent element, part)
    for (int i = 0; i < 4; i++)
      cm[i] = create_son_curv_map(e, i);
  }

  // create the four sons
  Element* sons[4];
  sons[0] = create_triangle(e->marker, e->vn[0], x0, x2, cm[0]);
  sons[1] = create_triangle(e->marker, x0, e->vn[1], x1, cm[1]);
  sons[2] = create_triangle(e->marker, x2, x1, e->vn[2], cm[2]);
  sons[3] = create_triangle(e->marker, x1, x2, x0, cm[3]);

  // update coefficients of curved reference mapping
  for (int i = 0; i < 4; i++)
    if (sons[i]->is_curved())
      sons[i]->cm->update_refmap_coefs(sons[i]);

  // deactivate this element and unregister from its nodes
  e->active = 0;
  nactive += 3;
  e->unref_all_nodes(this);
  // now the original edge nodes may no longer exist...

  // set correct boundary status and markers for the new nodes
  sons[0]->en[0]->bnd = bnd[0];  sons[0]->en[0]->marker = mrk[0];
  sons[0]->en[2]->bnd = bnd[2];  sons[0]->en[2]->marker = mrk[2];
  sons[1]->en[0]->bnd = bnd[0];  sons[1]->en[0]->marker = mrk[0];
  sons[1]->en[1]->bnd = bnd[1];  sons[1]->en[1]->marker = mrk[1];
  sons[2]->en[1]->bnd = bnd[1];  sons[2]->en[1]->marker = mrk[1];
  sons[2]->en[2]->bnd = bnd[2];  sons[2]->en[2]->marker = mrk[2];
  sons[3]->vn[0]->bnd = bnd[1];
  sons[3]->vn[1]->bnd = bnd[2];
  sons[3]->vn[2]->bnd = bnd[0];

  // copy son pointers (could not have been done earlier because of the union)
  memcpy(e->sons, sons, 4 * sizeof(Element*));
}


void Mesh::refine_quad(Element* e, int refinement)
{
  int i, j;
  Element* sons[4];

  // remember the markers of the edge nodes
  int bnd[4] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd,    e->en[3]->bnd    };
  int mrk[4] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker, e->en[3]->marker };

  // deactivate this element and unregister from its nodes
  e->active = false;
  nactive--;
  e->unref_all_nodes(this);
  // now the original edge nodes may no longer exist...

  CurvMap* cm[4];
  memset(cm, 0, sizeof(cm));

  // default refinement: one quad to four quads
  if (refinement == 0)
  {
    // obtain four mid-edge vertex nodes and one mid-element vetex node
    Node* x0 = get_vertex_node(e->vn[0]->id, e->vn[1]->id);
    Node* x1 = get_vertex_node(e->vn[1]->id, e->vn[2]->id);
    Node* x2 = get_vertex_node(e->vn[2]->id, e->vn[3]->id);
    Node* x3 = get_vertex_node(e->vn[3]->id, e->vn[0]->id);
    Node* mid = get_vertex_node(x0->id, x2->id);

    // adjust mid-edge coordinates if this is a curved element
    if (e->is_curved())
    {
      double2 pt[5] = { { 0.0,-1.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }, { -1.0, 0.0 }, { 0.0, 0.0 } };
      e->cm->get_mid_edge_points(e, pt, 5);
      x0->x = pt[0][0];  x0->y = pt[0][1];
      x1->x = pt[1][0];  x1->y = pt[1][1];
      x2->x = pt[2][0];  x2->y = pt[2][1];
      x3->x = pt[3][0];  x3->y = pt[3][1];
      mid->x = pt[4][0]; mid->y = pt[4][1];

      // create CurvMaps for sons (pointer to parent element, part)
      for (i = 0; i < 4; i++)
        cm[i] = create_son_curv_map(e, i);
    }

    // create the four sons
    sons[0] = create_quad(e->marker, e->vn[0], x0, mid, x3, cm[0]);
    sons[1] = create_quad(e->marker, x0, e->vn[1], x1, mid, cm[1]);
    sons[2] = create_quad(e->marker, mid, x1, e->vn[2], x2, cm[2]);
    sons[3] = create_quad(e->marker, x3, mid, x2, e->vn[3], cm[3]);
    nactive += 4;

    // set correct boundary markers for the new edge nodes
    for (i = 0; i < 4; i++)
    {
      j = (i > 0) ? i-1 : 3;
      sons[i]->en[j]->bnd = bnd[j];  sons[i]->en[j]->marker = mrk[j];
      sons[i]->en[i]->bnd = bnd[i];  sons[i]->en[i]->marker = mrk[i];
      sons[i]->vn[j]->bnd = bnd[j];
    }
  }
  // refinement '1': one quad to two 'horizontal' quads
  else if (refinement == 1)
  {
    Node* x1 = get_vertex_node(e->vn[1]->id, e->vn[2]->id);
    Node* x3 = get_vertex_node(e->vn[3]->id, e->vn[0]->id);

    // adjust mid-edge coordinates if this is a curved element
    if (e->is_curved())
    {
      double2 pt[2] = { { 1.0, 0.0 }, { -1.0, 0.0 } };
      e->cm->get_mid_edge_points(e, pt, 2);
      x1->x = pt[0][0];  x1->y = pt[0][1];
      x3->x = pt[1][0];  x3->y = pt[1][1];

      // create CurvMaps for sons (pointer to parent element, part)
      for (i = 0; i < 2; i++)
        cm[i] = create_son_curv_map(e, i + 4);
    }

    sons[0] = create_quad(e->marker, e->vn[0], e->vn[1], x1, x3, cm[0]);
    sons[1] = create_quad(e->marker, x3, x1, e->vn[2], e->vn[3], cm[1]);
    sons[2] = sons[3] = NULL;
    nactive += 2;

    sons[0]->en[0]->bnd = bnd[0];  sons[0]->en[0]->marker = mrk[0];
    sons[0]->en[1]->bnd = bnd[1];  sons[0]->en[1]->marker = mrk[1];
    sons[0]->en[3]->bnd = bnd[3];  sons[0]->en[3]->marker = mrk[3];
    sons[1]->en[1]->bnd = bnd[1];  sons[1]->en[1]->marker = mrk[1];
    sons[1]->en[2]->bnd = bnd[2];  sons[1]->en[2]->marker = mrk[2];
    sons[1]->en[3]->bnd = bnd[3];  sons[1]->en[3]->marker = mrk[3];
    sons[0]->vn[2]->bnd = bnd[1];
    sons[0]->vn[3]->bnd = bnd[3];
  }
  // refinement '2': one quad to two 'vertical' quads
  else if (refinement == 2)
  {
    Node* x0 = get_vertex_node(e->vn[0]->id, e->vn[1]->id);
    Node* x2 = get_vertex_node(e->vn[2]->id, e->vn[3]->id);

    // adjust mid-edge coordinates if this is a curved element
    if (e->is_curved())
    {
      double2 pt[2] = { { 0.0, -1.0 }, { 0.0, 1.0 } };
      e->cm->get_mid_edge_points(e, pt, 2);
      x0->x = pt[0][0];  x0->y = pt[0][1];
      x2->x = pt[1][0];  x2->y = pt[1][1];

      // create CurvMaps for sons (pointer to parent element, part)
      for (i = 0; i < 2; i++)
        cm[i] = create_son_curv_map(e, i + 6);
    }

    sons[0] = sons[1] = NULL;
    sons[2] = create_quad(e->marker, e->vn[0], x0, x2, e->vn[3], cm[0]);
    sons[3] = create_quad(e->marker, x0, e->vn[1], e->vn[2], x2, cm[1]);
    nactive += 2;

    sons[2]->en[0]->bnd = bnd[0];  sons[2]->en[0]->marker = mrk[0];
    sons[2]->en[2]->bnd = bnd[2];  sons[2]->en[2]->marker = mrk[2];
    sons[2]->en[3]->bnd = bnd[3];  sons[2]->en[3]->marker = mrk[3];
    sons[3]->en[0]->bnd = bnd[0];  sons[3]->en[0]->marker = mrk[0];
    sons[3]->en[1]->bnd = bnd[1];  sons[3]->en[1]->marker = mrk[1];
    sons[3]->en[2]->bnd = bnd[2];  sons[3]->en[2]->marker = mrk[2];
    sons[2]->vn[1]->bnd = bnd[0];
    sons[2]->vn[2]->bnd = bnd[2];
  }
  else assert(0);

  // update coefficients of curved reference mapping
  for (i = 0; i < 4; i++)
    if (sons[i] != NULL && sons[i]->cm != NULL)
      sons[i]->cm->update_refmap_coefs(sons[i]);

  // optimization: iro never gets worse
  if (e->iro_cache == 0)
    for (i = 0; i < 4; i++)
      if (sons[i] != NULL)
        sons[i]->iro_cache = 0;

  // copy son pointers (could not have been done earlier because of the union)
  memcpy(e->sons, sons, sizeof(sons));
}


void Mesh::unrefine_element_internal(Element* e)
{
  assert(!e->active);
  int i, s1, s2;

  // obtain markers and bnds from son elements
  int mrk[4], bnd[4];
  for (i = 0; i < e->nvert; i++)
  {
    get_edge_sons(e, i, s1, s2);
    assert(e->sons[s1]->active);
    mrk[i] = e->sons[s1]->en[i]->marker;
    bnd[i] = e->sons[s1]->en[i]->bnd;
  }

  // remove all sons
  for (i = 0; i < 4; i++)
  {
    Element* son = e->sons[i];
    if (son != NULL)
    {
      son->unref_all_nodes(this);
      if (son->cm != NULL) delete son->cm;
      elements.remove(son->id);
      nactive--;
    }
  }

  // recreate edge nodes
  for (i = 0; i < e->nvert; i++)
    e->en[i] = get_edge_node(e->vn[i]->id, e->vn[e->next_vert(i)]->id);

  e->ref_all_nodes();
  e->active = 1;
  nactive++;

  // restore edge node markers and bnds
  for (i = 0; i < e->nvert; i++)
  {
    e->en[i]->marker = mrk[i];
    e->en[i]->bnd = bnd[i];
  }
}


//// high-level element refinement /////////////////////////////////////////////////////////////////

void Mesh::refine_element(int id, int refinement)
{
  Element* e = get_element(id);
  if (!e->used) error("invalid element id number.");
  if (!e->active) error("attempt to refine element #%d which has been refined already.", e->id);

  if (e->is_triangle())
    refine_triangle(e);
  else
    refine_quad(e, refinement);

  seq = g_mesh_seq++;
}


void Mesh::refine_all_elements(int refinement)
{
  Element* e;
  elements.set_append_only(true);
  for_all_active_elements(e, this)
    refine_element(e->id, refinement);
  elements.set_append_only(false);
}


void Mesh::refine_by_criterion(int (*criterion)(Element*), int depth)
{
  Element* e;
  elements.set_append_only(true);
  for (int r, i = 0; i < depth; i++)
    for_all_active_elements(e, this)
      if ((r = criterion(e)) >= 0)
        refine_element(e->id, r);
  elements.set_append_only(false);
}


static int rtv_id;

static int rtv_criterion(Element* e)
{
  for (int i = 0; i < e->nvert; i++)
    if (e->vn[i]->id == rtv_id)
      return 0;
  return -1;
}

void Mesh::refine_towards_vertex(int vertex_id, int depth)
{
  rtv_id = vertex_id;
  refine_by_criterion(rtv_criterion, depth);
}


static int rtb_marker;
static bool rtb_aniso;
static char* rtb_vert;

static int rtb_criterion(Element* e)
{
  int i;
  for (i = 0; i < e->nvert; i++)
    if (e->en[i]->marker == rtb_marker || rtb_vert[e->vn[i]->id])
      break;

  if (i >= e->nvert) return -1;
  if (e->is_triangle() || !rtb_aniso) return 0;

  if ((e->en[0]->marker == rtb_marker && !rtb_vert[e->vn[2]->id] && !rtb_vert[e->vn[3]->id]) ||
      (e->en[2]->marker == rtb_marker && !rtb_vert[e->vn[0]->id] && !rtb_vert[e->vn[1]->id]) ||
      (e->en[0]->marker == rtb_marker && e->en[2]->marker == rtb_marker &&
       e->en[1]->marker != rtb_marker && e->en[3]->marker != rtb_marker)) return 1;

  if ((e->en[1]->marker == rtb_marker && !rtb_vert[e->vn[3]->id] && !rtb_vert[e->vn[0]->id]) ||
      (e->en[3]->marker == rtb_marker && !rtb_vert[e->vn[1]->id] && !rtb_vert[e->vn[2]->id]) ||
      (e->en[1]->marker == rtb_marker && e->en[3]->marker == rtb_marker &&
       e->en[0]->marker != rtb_marker && e->en[2]->marker != rtb_marker)) return 2;

  return 0;
}

void Mesh::refine_towards_boundary(int marker, int depth, bool aniso)
{
  rtb_marker = marker;
  rtb_aniso  = aniso;

  for (int i = 0; i < depth; i++)
  {
    int size = get_max_node_id()+1;
    rtb_vert = new char[size];
    memset(rtb_vert, 0, sizeof(char) * size);

    Element* e;
    for_all_active_elements(e, this)
      for (int j = 0; j < e->nvert; j++)
        if (e->en[j]->marker == marker)
          rtb_vert[e->vn[j]->id] = rtb_vert[e->vn[e->next_vert(j)]->id] = 1;

    refine_by_criterion(rtb_criterion, 1);
    delete [] rtb_vert;
  }
}


void Mesh::unrefine_element(int id)
{
  Element* e = get_element(id);
  if (!e->used) error("invalid element id number.");
  if (e->active) return;

  for (int i = 0; i < 4; i++)
    if (e->sons[i] != NULL)
      unrefine_element(e->sons[i]->id);

  unrefine_element_internal(e);
  seq = g_mesh_seq++;
}


void Mesh::unrefine_all_elements(bool keep_initial_refinements)
{
  // find inactive elements with active sons
  std::vector<int> list;
  Element* e;
  for_all_inactive_elements(e, this)
  {
    bool found = true;
    for (int i = 0; i < 4; i++)
      if (e->sons[i] != NULL && (!e->sons[i]->active ||
          (keep_initial_refinements && e->sons[i]->id < ninitial))  )
        { found = false; break; }

    if (found) list.push_back(e->id);
  }

  // unrefine the found elements
  for (int i = 0; i < list.size(); i++)
    unrefine_element(list[i]);
}
void Mesh::refine_triangle_to_quads(Element* e)
{
  // remember the markers of the edge nodes
  int bnd[3] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd    };
  int mrk[3] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker };

  // obtain three mid-edge and one gravity vertex nodes
  Node* x0 = get_vertex_node(e->vn[0]->id, e->vn[1]->id);
  Node* x1 = get_vertex_node(e->vn[1]->id, e->vn[2]->id);
  Node* x2 = get_vertex_node(e->vn[2]->id, e->vn[0]->id);
  Node* mid = get_vertex_node(x0->id, e->vn[1]->id);

  mid->x = (nodes[x0->id].x + nodes[x1->id].x + nodes[x2->id].x)/3;
  mid->y = (nodes[x0->id].y + nodes[x1->id].y + nodes[x2->id].y)/3;

  // adjust mid-edge and gravity coordinates if this is a curved element
  if (e->is_curved())
  {
    double2 pt[4] = { { 0.0,-1.0 }, { 0.0, 0.0 },{ -1.0, 0.0 }, { -0.33333333, -0.33333333 } };
    e->cm->get_mid_edge_points(e, pt, 4);
    x0->x = pt[0][0]; x0->y = pt[0][1];
    x1->x = pt[1][0]; x1->y = pt[1][1];
    x2->x = pt[2][0]; x2->y = pt[2][1];
    mid->x = pt[3][0]; mid->y = pt[3][1];
   }

  double angle2;
  int idx;
  CurvMap* cm[3];
  memset(cm, 0, sizeof(cm));

  // create CurvMaps for sons if this is a curved element
  if (e->is_curved())
  {
    for (idx = 0; idx < 2; idx++)
    {
      if (e->cm->nurbs[idx] != NULL)
      {
        cm[idx] = new CurvMap;
        memset(cm[idx], 0, sizeof(CurvMap));
        cm[idx+1] = new CurvMap;
        memset(cm[idx+1], 0, sizeof(CurvMap));
      }
    }

    idx=0;
    if (e->cm->nurbs[idx] != NULL)
    {
      angle2 = e->cm->nurbs[idx]->angle/2;
      Node* node_temp = get_vertex_node(e->vn[idx%3]->id, e->vn[(idx+1)%3]->id);

      for (int k = 0; k < 2; k++)
      {
        Node *en;
        int p1, p2;
        int idx2;

        if (k == 0)
        {
          p1 = e->vn[(idx)%3]->id;
          p2 = node_temp->id;
          if (idx == 0) idx2 = 0;
          if (idx == 1) idx2 = 1;
          if (idx == 2) continue;
        }
        else if (k == 1)
        {
          p1 = node_temp->id;
          p2 = e->vn[(idx+1)%3]->id;
          idx = (idx+1)%3;
          if (idx == 0) continue;
          if (idx == 1) idx2 = 0;
          if (idx == 2) idx2 = 0;
        }

        Nurbs* nurbs = new Nurbs;
        bool cricle = true;

        nurbs->arc = cricle;
        nurbs->degree = 2;

        int inner = 1, outer;
        inner = 1;
        nurbs->np = inner + 2;
        nurbs->pt = new double3[nurbs->np];

        nurbs->pt[0][0] = nodes[p1].x;
        nurbs->pt[0][1] = nodes[p1].y;
        nurbs->pt[0][2] = 1.0;

        nurbs->pt[inner+1][0] = nodes[p2].x;
        nurbs->pt[inner+1][1] = nodes[p2].y;
        nurbs->pt[inner+1][2] = 1.0;

        double angle = angle2;
        double a = (180.0 - angle) / 180.0 * M_PI;
        nurbs->angle = angle;

        // generate one control point
        double x = 1.0 / tan(a * 0.5);
        nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
        nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
        nurbs->pt[1][2] = cos((M_PI - a) * 0.5);

        int i;
        inner = 0;
        nurbs->nk = nurbs->degree + nurbs->np + 1;
        outer = nurbs->nk - inner;

        // knot vector is completed by 0.0 on the left and by 1.0 on the right
        nurbs->kv = new double[nurbs->nk];

        for (i = 0; i < outer/2; i++)
          nurbs->kv[i] = 0.0;
        for (i = outer/2 + inner; i < nurbs->nk; i++)
          nurbs->kv[i] = 1.0;
        nurbs->ref = 0;

        cm[idx]->toplevel = 1;
        cm[idx]->order = 4;
        cm[idx]->nurbs[idx2] = nurbs;
        nurbs->ref++;
      }
    }

    idx = 1;
    if (e->cm->nurbs[idx] != NULL)
    {
      angle2 = e->cm->nurbs[idx]->angle/2;
      Node* node_temp = get_vertex_node(e->vn[idx%3]->id, e->vn[(idx+1)%3]->id);
      for (int k = 0; k < 2; k++)
      {
        Node *en;
        int p1, p2;
        int idx2;
        if (k == 0)
        {
          p1 = e->vn[(idx)%3]->id;
          p2 = node_temp->id;
          if (idx == 0) idx2 = 0;
          if (idx == 1) idx2 = 1;
          if (idx == 2) continue;
        }
        else if (k == 1)
        {
          p1 = node_temp->id;
          p2 = e->vn[(idx+1)%3]->id;
          idx = (idx+1)%3;
          if (idx == 0) continue;
          if (idx == 1) idx2 = 0;
          if (idx == 2) idx2 = 0;
        }

        Nurbs* nurbs = new Nurbs;
        bool cricle = true;

        nurbs->arc = cricle;
        nurbs->degree = 2;
        int inner = 1, outer;
        inner = 1;
        nurbs->np = inner + 2;
        nurbs->pt = new double3[nurbs->np];

        nurbs->pt[0][0] = nodes[p1].x;
        nurbs->pt[0][1] = nodes[p1].y;
        nurbs->pt[0][2] = 1.0;

        nurbs->pt[inner+1][0] = nodes[p2].x;
        nurbs->pt[inner+1][1] = nodes[p2].y;
        nurbs->pt[inner+1][2] = 1.0;

        double angle = angle2;
        double a = (180.0 - angle) / 180.0 * M_PI;
        nurbs->angle = angle;

        // generate one control point
        double x = 1.0 / tan(a * 0.5);
        nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
        nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
        nurbs->pt[1][2] = cos((M_PI - a) * 0.5);

        int i;
        inner = 0;
        nurbs->nk = nurbs->degree + nurbs->np + 1;
        outer = nurbs->nk - inner;

        // knot vector is completed by 0.0 on the left and by 1.0 on the right
        nurbs->kv = new double[nurbs->nk];
        for (i = 0; i < outer/2; i++)
          nurbs->kv[i] = 0.0;
        for (i = outer/2 + inner; i < nurbs->nk; i++)
         nurbs->kv[i] = 1.0;
        nurbs->ref = 0;

        cm[idx]->toplevel = 1;
        cm[idx]->order = 4;
        cm[idx]->nurbs[idx2] = nurbs;
        nurbs->ref++;
      }
    }
  }

  // create the four sons
  Element* sons[4];
  sons[0] = create_quad(e->marker, e->vn[0], x0, mid, x2, cm[0]);
  sons[1] = create_quad(e->marker, x0, e->vn[1], x1, mid, cm[1]);
  sons[2] = create_quad(e->marker, x1, e->vn[2], x2, mid, cm[2]);
  sons[3] = NULL;

  // update coefficients of curved reference mapping
  for (int i = 0; i < 3; i++)
  {
    if (sons[i]->is_curved())
    {
      sons[i]->cm->update_refmap_coefs(sons[i]);
    }
  }

  // deactivate this element and unregister from its nodes
  e->active = 0;
  nactive += 3;
  e->unref_all_nodes(this);
  // now the original edge nodes may no longer exist...
  // set correct boundary status and markers for the new nodes
  sons[0]->en[0]->bnd = bnd[0];  sons[0]->en[0]->marker = mrk[0];
  sons[0]->en[3]->bnd = bnd[2];  sons[0]->en[3]->marker = mrk[2];
  sons[0]->vn[1]->bnd = bnd[0];

  sons[1]->en[0]->bnd = bnd[0];  sons[1]->en[0]->marker = mrk[0];
  sons[1]->en[1]->bnd = bnd[1];  sons[1]->en[1]->marker = mrk[1];
  sons[1]->vn[2]->bnd = bnd[1];

  sons[2]->en[0]->bnd = bnd[1];  sons[2]->en[0]->marker = mrk[1];
  sons[2]->en[1]->bnd = bnd[2];  sons[2]->en[1]->marker = mrk[2];
  sons[2]->vn[2]->bnd = bnd[2];

  // copy son pointers (could not have been done earlier because of the union)
  memcpy(e->sons, sons, 3 * sizeof(Element*));
}


void Mesh::refine_element_to_quads(int id, int refinement)
{
  Element* e = get_element(id);
  if (!e->used) error("invalid element id number.");
  if (!e->active) error("attempt to refine element #%d which has been refined already.", e->id);

  if (e->is_triangle())
    refine_triangle_to_quads(e);
  else
    refine_quad(e, refinement);

  seq = g_mesh_seq++;
}


void Mesh::refine_quad_to_triangles(Element* e)
{
  // remember the markers of the edge nodes
  int bnd[4] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd,    e->en[3]->bnd };
  int mrk[4] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker, e->en[3]->marker };


  bool bcheck = true;  ///< if bcheck is true, it is default add a new edge between
                       ///<  vn[0] and vn[2]
  double length_x_0_2 = (e->vn[0]->x - e->vn[2]->x)*(e->vn[0]->x - e->vn[2]->x);
  double length_x_1_3 = (e->vn[1]->x - e->vn[3]->x)*(e->vn[1]->x - e->vn[3]->x);

  double length_y_0_2 = (e->vn[0]->y - e->vn[2]->y)*(e->vn[0]->y - e->vn[2]->y);
  double length_y_1_3 = (e->vn[1]->y - e->vn[3]->y)*(e->vn[1]->y - e->vn[3]->y);

  if ((length_x_0_2 + length_y_0_2) > (length_x_1_3 + length_y_1_3))
  {
    bcheck = false;
  }

  double angle2;
  int idx;
  CurvMap* cm[2];
  memset(cm, 0, sizeof(cm));

  // create CurvMaps for sons if this is a curved element
  if (e->is_curved())
  {
    int i_case2 = 0;
    if (bcheck == true)
    {
      if ((e->cm->nurbs[0] != NULL) || (e->cm->nurbs[1] != NULL))
      {
        cm[0] = new CurvMap;
        memset(cm[0], 0, sizeof(CurvMap));
      }
      if ((e->cm->nurbs[2] != NULL) || (e->cm->nurbs[3] != NULL))
      {
        cm[1] = new CurvMap;
        memset(cm[1], 0, sizeof(CurvMap));
      }
    }
    else if (bcheck == false)
    {
      if ((e->cm->nurbs[1] != NULL) || (e->cm->nurbs[2] != NULL))
      {
        cm[0] = new CurvMap;
        memset(cm[0], 0, sizeof(CurvMap));
      }
      if ((e->cm->nurbs[3] != NULL) || (e->cm->nurbs[0] != NULL))
      {
        cm[1] = new CurvMap;
        memset(cm[1], 0, sizeof(CurvMap));
      }
      i_case2 = 1; //switch to the shorter diagonal
    }

    for (int k = 0; k < 2; k++)
    {
      for (idx = 0 + 2*k; idx < 2 + 2*k; idx++)
      {
        if (e->cm->nurbs[(idx + i_case2)%4] != NULL)
        {
          angle2 = e->cm->nurbs[(idx + i_case2)%4]->angle;

          Node *en;
          int p1, p2;
          int idx2 = idx;

          p1 = e->vn[(idx + i_case2)%4]->id;
          p2 = e->vn[(idx + i_case2 + 1)%4]->id;  //node_temp->id;

          Nurbs* nurbs = new Nurbs;
          bool cricle = true;

          nurbs->arc = cricle;
          nurbs->degree = 2;

          int inner = 1, outer;
          inner = 1;
          nurbs->np = inner + 2;
          nurbs->pt = new double3[nurbs->np];

          nurbs->pt[0][0] = nodes[p1].x;
          nurbs->pt[0][1] = nodes[p1].y;
          nurbs->pt[0][2] = 1.0;

          nurbs->pt[inner+1][0] = nodes[p2].x;
          nurbs->pt[inner+1][1] = nodes[p2].y;
          nurbs->pt[inner+1][2] = 1.0;

          double angle = angle2;
          double a = (180.0 - angle) / 180.0 * M_PI;
          nurbs->angle = angle;

          // generate one control point
          double x = 1.0 / tan(a * 0.5);
          nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
          nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
          nurbs->pt[1][2] = cos((M_PI - a) * 0.5);

          int i;
          inner = 0;
          nurbs->nk = nurbs->degree + nurbs->np + 1;
          outer = nurbs->nk - inner;

          // knot vector is completed by 0.0 on the left and by 1.0 on the right
          nurbs->kv = new double[nurbs->nk];

          for (i = 0; i < outer/2; i++)
            nurbs->kv[i] = 0.0;
          for (i = outer/2 + inner; i < nurbs->nk; i++)
            nurbs->kv[i] = 1.0;
          nurbs->ref = 0;

          cm[k]->toplevel = 1;
          cm[k]->order = 4;
          cm[k]->nurbs[idx%2] = nurbs;
          nurbs->ref++;
        }
      }
    }
  }

  // create the four sons
  Element* sons[4];
  if (bcheck == true)
  {
    sons[0] = create_triangle(e->marker, e->vn[0], e->vn[1], e->vn[2], cm[0]);
    sons[1] = create_triangle(e->marker, e->vn[2], e->vn[3], e->vn[0], cm[1]);
    sons[2] = NULL; //create_quad(e->marker, x1, e->vn[2], x2, mid, cm[2]);
    sons[3] = NULL;
  }
  else
  {
    sons[0] = create_triangle(e->marker, e->vn[1], e->vn[2], e->vn[3], cm[0]);
    sons[1] = create_triangle(e->marker, e->vn[3], e->vn[0], e->vn[1], cm[1]);
    sons[2] = NULL; //create_quad(e->marker, x1, e->vn[2], x2, mid, cm[2]);
    sons[3] = NULL;
  }

  // update coefficients of curved reference mapping
  for (int i = 0; i < 2; i++)
  {
    if (sons[i]->is_curved())
    {
      sons[i]->cm->update_refmap_coefs(sons[i]);
    }
  }

  // deactivate this element and unregister from its nodes
  e->active = 0;
  nactive += 3;
  e->unref_all_nodes(this);
  // now the original edge nodes may no longer exist...
  // set correct boundary status and markers for the new nodes
  if (bcheck == true)
  {
    sons[0]->en[0]->bnd = bnd[0];  sons[0]->en[0]->marker = mrk[0];
    sons[0]->en[1]->bnd = bnd[1];  sons[0]->en[1]->marker = mrk[1];
    sons[0]->vn[1]->bnd = bnd[0];

    sons[1]->en[0]->bnd = bnd[2];  sons[1]->en[0]->marker = mrk[2];
    sons[1]->en[1]->bnd = bnd[3];  sons[1]->en[1]->marker = mrk[3];
    sons[1]->vn[2]->bnd = bnd[1];
  }
  else
  {
    sons[0]->en[0]->bnd = bnd[1];  sons[0]->en[0]->marker = mrk[1];
    sons[0]->en[1]->bnd = bnd[2];  sons[0]->en[1]->marker = mrk[2];
    sons[0]->vn[1]->bnd = bnd[1];

    sons[1]->en[0]->bnd = bnd[3];  sons[1]->en[0]->marker = mrk[3];
    sons[1]->en[1]->bnd = bnd[0];  sons[1]->en[1]->marker = mrk[0];
    sons[1]->vn[2]->bnd = bnd[0];
  }

  // copy son pointers (could not have been done earlier because of the union)
  memcpy(e->sons, sons, 2 * sizeof(Element*));
}


void Mesh::refine_element_to_triangles(int id)
{
  Element* e = get_element(id);
  if (!e->used) error("invalid element id number.");
  if (!e->active) error("attempt to refine element #%d which has been refined already.", e->id);

  if (e->is_triangle())
    refine_triangle(e);
  else
    refine_quad_to_triangles(e);

  seq = g_mesh_seq++;
}
