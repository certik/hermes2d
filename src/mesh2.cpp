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
#include "mesh_parser.h"


extern unsigned g_mesh_seq;


//// load //////////////////////////////////////////////////////////////////////////////////////////

/*************** OLD OLD OLD ***************************************************************/

#define eof_error error("premature end of file")

char* Mesh::get_line(FILE* f)
{
  static char line[1000];

  // read one line, skipping empty ones and those starting with '*' or '#'
  while (fgets(line, 1000-1, f) != NULL)
  {
    char* p = line;
    while (*p && (unsigned) *p <= ' ') p++;

    if (*p && *p != '*' && *p != '#') return line;
  }

  return NULL;
}


Nurbs* Mesh::load_nurbs(FILE* f, Node** en, int &p1, int &p2)
{
  int i;
  char* line;
  Nurbs* nurbs = new Nurbs;

  // read the end point indices
  if ((line = get_line(f)) == NULL) eof_error;
  if (sscanf(line, "%d %d", &p1, &p2) != 2) error("error reading curved boundary data (end point indices)");
  *en = peek_edge_node(p1, p2);
  if (*en == NULL) error("error reading curved boundary data (edge %d-%d does not exist)", p1, p2);

  // degree of curved edge
  if ((line = get_line(f)) == NULL) eof_error;
  if (sscanf(line, "%d", &(nurbs->degree)) != 1)
    error("error reading curved boundary data for edge %d-%d (degree)", p1, p2);

  // create a circular arc if degree == 0
  bool circle = (nurbs->degree == 0);
  if (circle) nurbs->degree = 2;
  nurbs->arc = circle;

  // load control points of curved edge
  int inner = 1, outer;
  if (!circle)
  {
    if ((line = get_line(f)) == NULL) eof_error;
    if (sscanf(line, "%d", &inner) != 1)
      error("error reading curved boundary data for edge %d-%d (# of control points)", p1, p2);
  }
  nurbs->np = inner + 2;

  // edge endpoints are also control points, with weight 1.0
  nurbs->pt = new double3[nurbs->np];
  nurbs->pt[0][0] = nodes[p1].x;
  nurbs->pt[0][1] = nodes[p1].y;
  nurbs->pt[0][2] = 1.0;
  nurbs->pt[inner+1][0] = nodes[p2].x;
  nurbs->pt[inner+1][1] = nodes[p2].y;
  nurbs->pt[inner+1][2] = 1.0;

  if (!circle)
  {
    // load inner control points
    for (i = 1; i <= inner; i++)
    {
      if ((line = get_line(f)) == NULL) eof_error;
      if (sscanf(line, "%lf %lf %lf", &(nurbs->pt[i][0]), &(nurbs->pt[i][1]), &(nurbs->pt[i][2])) != 3)
        error("error reading curved boundary data for edge %d-%d (control points)", p1, p2);
    }
  }
  else
  {
    // load the arc angle
    if ((line = get_line(f)) == NULL) eof_error;
    double angle;
    if (sscanf(line, "%lf", &angle) != 1)
      error("error reading curved boundary data for edge %d-%d (arc angle)", p1, p2);
    nurbs->angle = angle;
    angle = (180.0 - angle) / 180.0 * M_PI;

    // generate one control point
    double x = 1.0 / tan(angle * 0.5);
    nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
    nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
    nurbs->pt[1][2] = cos((M_PI - angle) * 0.5);
  }

  // load the number of knot vector points
  inner = 0;
  if (!circle)
  {
    if ((line = get_line(f)) == NULL) eof_error;
    if (sscanf(line, "%d", &(inner)) != 1)
      error("error reading curved boundary data for edge %d-%d (# of knot vector points)", p1, p2);
  }

  nurbs->nk = nurbs->degree + nurbs->np + 1;
  outer = nurbs->nk - inner;
  if (outer & 1 == 1)
    error("error reading curved boundary data for edge %d-%d (wrong number of knot points)", p1, p2);

  // knot vector is completed by 0.0 on the left and by 1.0 on the right
  nurbs->kv = new double[nurbs->nk];
  for (i = 0; i < outer/2; i++)
    nurbs->kv[i] = 0.0;
  for (i = outer/2; i < inner + outer/2; i++)
  {
    if ((line = get_line(f)) == NULL) eof_error;
    if (sscanf(line, "%lf", &(nurbs->kv[i])) != 1)
      error("error reading curved boundary data for edge %d-%d (knot points)", p1, p2);
  }
  for (i = outer/2 + inner; i < nurbs->nk; i++)
    nurbs->kv[i] = 1.0;

  nurbs->ref = 0;
  return nurbs;
}


/// Returns a NURBS curve with reversed control points and inverted knot vector.
/// Used for curved edges inside a mesh, where two mirror Nurbs have to be created
/// for the adjacent elements
///
Nurbs* Mesh::reverse_nurbs(Nurbs* nurbs)
{
  Nurbs* rev = new Nurbs;
  *rev = *nurbs;
  rev->twin = true;

  rev->pt = new double3[nurbs->np];
  for (int i = 0; i < nurbs->np; i++)
  {
    rev->pt[nurbs->np-1 - i][0] = nurbs->pt[i][0];
    rev->pt[nurbs->np-1 - i][1] = nurbs->pt[i][1];
    rev->pt[nurbs->np-1 - i][2] = nurbs->pt[i][2];
  }

  rev->kv = new double[nurbs->nk];
  for (int i = 0; i < nurbs->nk; i++)
    rev->kv[i] = nurbs->kv[i];
  for (int i = nurbs->degree + 1; i < nurbs->nk - nurbs->degree - 1; i++)
    rev->kv[nurbs->nk-1 - i] = 1.0 - nurbs->kv[i];

  rev->arc = nurbs->arc;
  rev->angle = -nurbs->angle;
  return rev;
}

// computing vector length
double vector_length(double a_1, double a_2)
{
  return sqrt(sqr(a_1) + sqr(a_2));
}

// checking whether the points p, q, r lie on the same line
bool same_line(double p_1, double p_2, double q_1, double q_2, double r_1, double r_2)
{
  double pq_1 = q_1 - p_1, pq_2 = q_2 - p_2, pr_1 = r_1 - p_1, pr_2 = r_2 - p_2;
  double length_pq = vector_length(pq_1, pq_2);
  double length_pr = vector_length(pr_1, pr_2);
  double sin_angle = (pq_1*pr_2 - pq_2*pr_1)/(length_pq*length_pr);
  if(fabs(sin_angle) < 1e-8) return true;
  else return false;
}

// checking whether the angle of vectors 'a' and 'b' is between zero and Pi
bool is_convex(double a_1, double a_2, double b_1, double b_2)
{
  if(a_1*b_2 - a_2*b_1 > 0) return true;
  else return false;
}

void check_triangle(int i, Node *&v0, Node *&v1, Node *&v2)
{
  // checking that all edges have nonzero length
  double
    length_1 = vector_length(v1->x - v0->x, v1->y - v0->y),
    length_2 = vector_length(v2->x - v1->x, v2->y - v1->y),
    length_3 = vector_length(v0->x - v2->x, v0->y - v2->y);
  if(length_1 < 1e-14 || length_2 < 1e-14 || length_3 < 1e-14)
    error("Edge of triangular element #%d has length less than 1e-14.", i);

  // checking that vertices do not lie on the same line
  if(same_line(v0->x, v0->y, v1->x, v1->y, v2->x, v2->y))
    error("Triangular element #%d: all vertices lie on the same line.", i);

  // checking positive orientation. If not positive, swapping vertices
  if (!is_convex(v1->x - v0->x, v1->y - v0->y, v2->x - v0->x, v2->y - v0->y)) {
    warn("Triangular element #%d not positively oriented, swapping vertices.", i);
    std::swap(v1, v2);
  }
}

void check_quad(int i, Node *&v0, Node *&v1, Node *&v2, Node *&v3)
{
  // checking that all edges have nonzero length
  double
    length_1 = vector_length(v1->x - v0->x, v1->y - v0->y),
    length_2 = vector_length(v2->x - v1->x, v2->y - v1->y),
    length_3 = vector_length(v3->x - v2->x, v3->y - v2->y),
    length_4 = vector_length(v0->x - v3->x, v0->y - v3->y);
  if(length_1 < 1e-14 || length_2 < 1e-14 || length_3 < 1e-14 || length_4 < 1e-14)
    error("Edge of quad element #%d has length less than 1e-14.", i);

  // checking that both diagonals have nonzero length
  double
    diag_1 = vector_length(v2->x - v0->x, v2->y - v0->y),
    diag_2 = vector_length(v3->x - v1->x, v3->y - v1->y);
  if(diag_1 < 1e-14 || diag_2 < 1e-14)
    error("Diagonal of quad element #%d has length less than 1e-14.", i);

  // checking that vertices v0, v1, v2 do not lie on the same line
  if(same_line(v0->x, v0->y, v1->x, v1->y, v2->x, v2->y))
    error("Quad element #%d: vertices v0, v1, v2 lie on the same line.", i);
  // checking that vertices v0, v1, v3 do not lie on the same line
  if(same_line(v0->x, v0->y, v1->x, v1->y, v3->x, v3->y))
    error("Quad element #%d: vertices v0, v1, v3 lie on the same line.", i);
  // checking that vertices v0, v2, v3 do not lie on the same line
  if(same_line(v0->x, v0->y, v2->x, v2->y, v3->x, v3->y))
    error("Quad element #%d: vertices v0, v2, v3 lie on the same line.", i);
  // checking that vertices v1, v2, v3 do not lie on the same line
  if(same_line(v1->x, v1->y, v2->x, v2->y, v3->x, v3->y))
    error("Quad element #%d: vertices v1, v2, v3 lie on the same line.", i);

  // checking that vertex v1 lies on the right of the diagonal v2-v0
  int vertex_1_ok = is_convex(v1->x - v0->x, v1->y - v0->y, v2->x - v0->x, v2->y - v0->y);
  if(!vertex_1_ok) error("Vertex v1 of quad element #%d does not lie on the right of the diagonal v2-v0.", i);
  // checking that vertex v3 lies on the left of the diagonal v2-v0
  int vertex_3_ok = is_convex(v2->x - v0->x, v2->y - v0->y, v3->x - v0->x, v3->y - v0->y);
  if(!vertex_3_ok) error("Vertex v3 of quad element #%d does not lie on the left of the diagonal v2-v0.", i);
  // checking that vertex v2 lies on the right of the diagonal v3-v1
  int vertex_2_ok = is_convex(v2->x - v1->x, v2->y - v1->y, v3->x - v1->x, v3->y - v1->y);
  if(!vertex_2_ok) error("Vertex v2 of quad element #%d does not lie on the right of the diagonal v3-v1.", i);
  // checking that vertex v0 lies on the left of the diagonal v3-v1
  int vertex_0_ok = is_convex(v3->x - v1->x, v3->y - v1->y, v0->x - v1->x, v0->y - v1->y);
  if(!vertex_0_ok) error("Vertex v0 of quad element #%d does not lie on the left of the diagonal v2-v1.", i);
}

void Mesh::load_old(const char* filename)
{
  // open the mesh file
  FILE* f = fopen(filename, "r");
  if (f == NULL) error("could not open the mesh file %s", filename);
  this->load_stream(f);
}

void Mesh::load_str(char* mesh)
{
  // open the mesh file
  FILE* f = fmemopen(mesh, strlen(mesh), "r");
  if (f == NULL) error("could not create the read buffer");
  this->load_stream(f);
}

/*
   Loads the mesh from a stream.
*/
void Mesh::load_stream(FILE *f)
{
  int i, j, k, n, maj, min;
  char* line;

  // check file version
  if ((line = get_line(f)) == NULL) eof_error;
  if (sscanf(line, "%d %d", &maj, &min) != 2) error("could not read file version");
  if (maj > 1) error("unsupported file version");

  // read the number of vertices
  if ((line = get_line(f)) == NULL) eof_error;
  if (sscanf(line, "%d", &n) != 1) error("could not read the number of vertices");

  // free all current data
  free();

  // create a hash table large enough
  int size = DEFAULT_HASH_SIZE;
  while (size < 8*n) size *= 2;
  HashTable::init(size);

  // load vertices: create top-level vertex nodes
  for (i = 0; i < n; i++)
  {
    Node* node = nodes.add();
    assert(node->id == i);
    node->ref = TOP_LEVEL_REF;
    node->type = TYPE_VERTEX;
    node->bnd = 0;
    node->p1 = node->p2 = -1;
    node->next_hash = NULL;

    if ((line = get_line(f)) == NULL) eof_error;
    if (sscanf(line, "%lf %lf", &node->x, &node->y) != 2) error("error reading vertex data");
  }
  ntopvert = n;

  // read the number of elements
  if ((line = get_line(f)) == NULL) eof_error;
  if (sscanf(line, "%d", &n) != 1) error("could not read the number of elements");

  // load elements
  for (i = 0; i < n; i++)
  {
    if ((line = get_line(f)) == NULL) eof_error;

    int ret, idx[5];
    if ((ret = sscanf(line, "%d %d %d %d %d", idx, idx+1, idx+2, idx+3, idx+4)) != 4 && ret != 5)
      error("error reading elements");

    for (j = 0; j < ret-1; j++)
      if (idx[j] < 0 || idx[j] >= ntopvert)
        error("error reading elements: node %d does not exist", idx[j]);

    Node *v0 = &nodes[idx[0]], *v1 = &nodes[idx[1]], *v2 = &nodes[idx[2]];
    if (ret == 4)
    {
      check_triangle(i, v0, v1, v2);
      create_triangle(idx[3], v0, v1, v2, NULL);
    }
    else
    {
      Node *v3 = &nodes[idx[3]];
      check_quad(i, v0, v1, v2, v3);
      create_quad(idx[4], v0, v1, v2, v3, NULL);
    }
  }
  nbase = nactive = n;

  // read the number of boundary data
  if ((line = get_line(f)) == NULL) eof_error;
  if (sscanf(line, "%d", &n) != 1) error("could not read the number of boundary markers\n");

  // load boundary data
  Node* en;
  for (i = 0; i < n; i++)
  {
    if ((line = get_line(f)) == NULL) eof_error;

    int v1, v2, marker;
    if (sscanf(line, "%d %d %d", &v1, &v2, &marker) != 3) error("error reading boundary marker data");

    en = peek_edge_node(v1, v2);
    if (en == NULL) error("boundary data error (edge %d-%d does not exist)", v1, v2);
    en->marker = marker;

    if (marker > 0)
    {
      nodes[v1].bnd = 1;
      nodes[v2].bnd = 1;
      en->bnd = 1;
    }
  }

  // check that all boundary edges have a marker assigned
  for_all_edge_nodes(en, this)
    if (en->ref < 2 && en->marker == 0)
      warn("boundary edge node does not have a boundary marker");

  // read the number of curved edges
  if ((line = get_line(f)) == NULL) eof_error;
  if (sscanf(line, "%d", &n) != 1) error("could not read the number of curved edges");

  // load curved edges
  for (i = 0; i < n; i++)
  {
    // load the control points, knot vector, etc.
    Node* en;
    int p1, p2;
    Nurbs* nurbs = load_nurbs(f, &en, p1, p2);

    // assign the nurbs to the elements sharing the edge node
    for (k = 0; k < 2; k++)
    {
      Element* e = en->elem[k];
      if (e == NULL) continue;

      if (e->cm == NULL)
      {
        e->cm = new CurvMap;
        memset(e->cm, 0, sizeof(CurvMap));
        e->cm->toplevel = 1;
        e->cm->order = 4;
      }

      int idx = -1;
      for (j = 0; j < e->nvert; j++)
        if (e->en[j] == en) { idx = j; break; }
      assert(idx >= 0);

      if (e->vn[idx]->id == p1)
      {
        e->cm->nurbs[idx] = nurbs;
        nurbs->ref++;
      }
      else
      {
        Nurbs* nurbs_rev = reverse_nurbs(nurbs);
        e->cm->nurbs[idx] = nurbs_rev;
        nurbs_rev->ref++;
      }
    }
    if (!nurbs->ref) delete nurbs;
  }

  // read the number of initial refinements
  //if () eof_error;
  //if (sscanf(line, "%d", &n) != 1) error("could not read the number of initial refinements");

  if ((line = get_line(f)) == NULL ||
      sscanf(line, "%d", &n) != 1)
  {
    warn("could not read the number of initial refinements");
  }
  else
  {
    // perform initial refinements
    for (i = 0; i < n; i++)
    {
      if ((line = get_line(f)) == NULL) eof_error;
      int id, ref;
      if (sscanf(line, "%d %d", &id, &ref) != 2)
        error("error reading initial refinement data");
      refine_element(id, ref);
    }
  }
  ninitial = elements.get_num_items();

  // update refmap coefs of curvilinear elements
  Element* e;
  for_all_elements(e, this)
    if (e->cm != NULL)
      e->cm->update_refmap_coefs(e);

  fclose(f);
  seq = g_mesh_seq++;
}

/*************** OLD OLD OLD ***************************************************************/


//// load_new //////////////////////////////////////////////////////////////////////////////////////

void Mesh::load(const char* filename, bool debug)
{
  int i, j, k, n;
  Node* en;

  // open the mesh file
  FILE* f = fopen(filename, "r");
  if (f == NULL) error("could not open the mesh file %s", filename);

  // free all current data
  free();

  // parse the file
  mesh_parser_init(f, filename);
  mesh_parser_run(debug);
  fclose(f);

  //// vertices ////////////////////////////////////////////////////////////////

  MSymbol* sym = mesh_parser_find_symbol("vertices");
  if (sym == NULL) error("%s: 'vertices' not found.", filename);
  n = sym->data->n;
  if (n < 0) error("%s: 'vertices' must be a list.", filename);
  if (n < 2) error("%s: invalid number of vertices.", filename);

  // create a hash table large enough
  int size = DEFAULT_HASH_SIZE;
  while (size < 8*n) size *= 2;
  HashTable::init(size);

  // create top-level vertex nodes
  MItem* pair = sym->data->list;
  for (i = 0; i < n; i++, pair = pair->next)
  {
    Node* node = nodes.add();
    assert(node->id == i);
    node->ref = TOP_LEVEL_REF;
    node->type = TYPE_VERTEX;
    node->bnd = 0;
    node->p1 = node->p2 = -1;
    node->next_hash = NULL;

    if (!mesh_parser_get_doubles(pair, 2, &node->x, &node->y))
      error("%s: invalid vertex #%d.", filename, i);
  }
  ntopvert = n;

  //// elements ////////////////////////////////////////////////////////////////

  sym = mesh_parser_find_symbol("elements");
  if (sym == NULL) error("%s: 'elements' not found.", filename);
  n = sym->data->n;
  if (n < 0) error("%s: 'elements' must be a list.", filename);
  if (n < 1) error("%s: no elements defined.", filename);

  // create elements
  MItem* elem = sym->data->list;
  nactive = 0;
  for (i = 0; i < n; i++, elem = elem->next)
  {
    // read and check vertex indices
    int nv = elem->n, idx[5];
    if (!nv) { elements.skip_slot(); continue; }
    if (nv < 4 || nv > 5)
      error("%s: element #%d: wrong number of vertex indices.", filename, i);
    if (!mesh_parser_get_ints(elem, nv, &idx[0], &idx[1], &idx[2], &idx[3], &idx[4]))
      error("%s: invalid definition of element #%d.", filename, i);
    for (j = 0; j < nv-1; j++)
      if (idx[j] < 0 || idx[j] >= ntopvert)
        error("%s: error creating element #%d: vertex #%d does not exist.", filename, i, idx[j]);

    // create triangle/quad
    Node *v0 = &nodes[idx[0]], *v1 = &nodes[idx[1]], *v2 = &nodes[idx[2]];
    if (nv == 4)
    {
      check_triangle(i, v0, v1, v2);
      create_triangle(idx[3], v0, v1, v2, NULL);
    }
    else
    {
      Node *v3 = &nodes[idx[3]];
      check_quad(i, v0, v1, v2, v3);
      create_quad(idx[4], v0, v1, v2, v3, NULL);
    }
    nactive++;
  }
  nbase = n;

  //// boundaries //////////////////////////////////////////////////////////////

  sym = mesh_parser_find_symbol("boundaries");
  if (sym != NULL)
  {
    n = sym->data->n;
    if (n < 0) error("%s: 'boundaries' must be a list.", filename);

    // read boundary data
    MItem* triple = sym->data->list;
    for (i = 0; i < n; i++, triple = triple->next)
    {
      int v1, v2, marker;
      if (!mesh_parser_get_ints(triple, 3, &v1, &v2, &marker))
        error("%s: invalid boundary data #%d.", filename, i);

      en = peek_edge_node(v1, v2);
      if (en == NULL)
        error("%s: boundary data #%d: edge %d-%d does not exist", filename, i, v1, v2);
      en->marker = marker;

      if (marker > 0)
      {
        nodes[v1].bnd = 1;
        nodes[v2].bnd = 1;
        en->bnd = 1;
      }
    }
  }

  // check that all boundary edges have a marker assigned
  for_all_edge_nodes(en, this)
    if (en->ref < 2 && en->marker == 0)
      warn("boundary edge node does not have a boundary marker");

  //// curves //////////////////////////////////////////////////////////////////

  sym = mesh_parser_find_symbol("curves");
  if (sym != NULL)
  {
    n = sym->data->n;
    if (n < 0) error("%s: 'curves' must be a list.", filename);

    // load curved edges
    MItem* curve = sym->data->list;
    for (i = 0; i < n; i++, curve = curve->next)
    {
      // load the control points, knot vector, etc.
      Node* en;
      int p1, p2;
      Nurbs* nurbs = load_nurbs_new(curve, i, &en, p1, p2);

      // assign the nurbs to the elements sharing the edge node
      for (k = 0; k < 2; k++)
      {
        Element* e = en->elem[k];
        if (e == NULL) continue;

        if (e->cm == NULL)
        {
          e->cm = new CurvMap;
          memset(e->cm, 0, sizeof(CurvMap));
          e->cm->toplevel = 1;
          e->cm->order = 4;
        }

        int idx = -1;
        for (j = 0; j < e->nvert; j++)
          if (e->en[j] == en) { idx = j; break; }
        assert(idx >= 0);

        if (e->vn[idx]->id == p1)
        {
          e->cm->nurbs[idx] = nurbs;
          nurbs->ref++;
        }
        else
        {
          Nurbs* nurbs_rev = reverse_nurbs(nurbs);
          e->cm->nurbs[idx] = nurbs_rev;
          nurbs_rev->ref++;
        }
      }
      if (!nurbs->ref) delete nurbs;
    }
  }

  // update refmap coefs of curvilinear elements
  Element* e;
  for_all_elements(e, this)
    if (e->cm != NULL)
      e->cm->update_refmap_coefs(e);

  //// refinements /////////////////////////////////////////////////////////////

  sym = mesh_parser_find_symbol("refinements");
  if (sym != NULL)
  {
    n = sym->data->n;
    if (n < 0) error("%s: 'refinements' must be a list.", filename);

    // perform initial refinements
    MItem* pair = sym->data->list;
    for (i = 0; i < n; i++, pair = pair->next)
    {
      int id, ref;
      if (!mesh_parser_get_ints(pair, 2, &id, &ref))
        error("%s: invalid refinement #%d.", filename, i);
      refine_element(id, ref);
    }
  }
  ninitial = elements.get_num_items();

  mesh_parser_free();
  seq = g_mesh_seq++;
}


//// load_nurbs ////////////////////////////////////////////////////////////////////////////////////

Nurbs* Mesh::load_nurbs_new(MItem* curve, int id, Node** en, int &p1, int &p2)
{
  int i;
  Nurbs* nurbs = new Nurbs;

  if (curve == NULL || curve->n < 0 || curve->n != 3 && curve->n != 5)
    error("invalid curve #%d.", id);
  bool circle = (curve->n == 3);
  nurbs->arc = circle;

  // read the end point indices
  MItem* edge = curve->list;
  if (edge->n >= 0 || !is_int(edge->val))
    error("curve #%d: invalid edge definition.", id);
  p1 = (int) edge->val;
  edge = edge->next;

  if (edge->n >= 0 || !is_int(edge->val))
    error("curve #%d: invalid edge definition.", id);
  p2 = (int) edge->val;
  edge = edge->next;

  *en = peek_edge_node(p1, p2);
  if (*en == NULL)
    error("curve #%d: edge %d-%d does not exist.", id, p1, p2);

  // degree of curved edge
  MItem* deg = edge;
  nurbs->degree = 2;
  if (!circle)
  {
    if (deg == NULL || deg->n >= 0 || !is_int(deg->val) || deg->val < 0 || deg->val == 1)
      error("curve #%d: invalid degee.", id);
    nurbs->degree = (int) deg->val;
  }

  // get the number of control points
  MItem* pts = deg->next;
  int inner = 1, outer;
  if (!circle)
  {
    if (pts == NULL || pts->n < 0)
      error("curve #%d: control points not defined.", id);
    inner = pts->n;
  }
  nurbs->np = inner + 2;

  // edge endpoints are also control points, with weight 1.0
  nurbs->pt = new double3[nurbs->np];
  nurbs->pt[0][0] = nodes[p1].x;
  nurbs->pt[0][1] = nodes[p1].y;
  nurbs->pt[0][2] = 1.0;
  nurbs->pt[inner+1][0] = nodes[p2].x;
  nurbs->pt[inner+1][1] = nodes[p2].y;
  nurbs->pt[inner+1][2] = 1.0;

  if (!circle)
  {
    // read inner control points
    MItem* it = pts->list;
    for (i = 1; i <= inner; i++, it = it->next)
    {
      if (!mesh_parser_get_doubles(it, 3, &(nurbs->pt[i][0]), &(nurbs->pt[i][1]), &(nurbs->pt[i][2])))
        error("curve #%d: invalid control point #%d.", id, i-1);
    }
  }
  else
  {
    // read the arc angle
    MItem* angle = deg;
    if (angle == NULL || angle->n >= 0)
      error("curve #%d: invalid arc angle.", id);
    double a = (180.0 - angle->val) / 180.0 * M_PI;
    nurbs->angle = angle->val;

    // generate one control point
    double x = 1.0 / tan(a * 0.5);
    nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
    nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
    nurbs->pt[1][2] = cos((M_PI - a) * 0.5);
  }

  // get the number of knot vector points
  inner = 0;
  MItem* knot;
  if (!circle && knot != NULL)
  {
    knot = pts->next;
    if (knot->n < 0) error("curve #%d: invalid knot vector.", id);
    inner = knot->n;
  }

  nurbs->nk = nurbs->degree + nurbs->np + 1;
  outer = nurbs->nk - inner;
  if (outer & 1 == 1)
    error("curve #%d: incorrect number of knot points.", id);

  // knot vector is completed by 0.0 on the left and by 1.0 on the right
  nurbs->kv = new double[nurbs->nk];
  for (i = 0; i < outer/2; i++)
    nurbs->kv[i] = 0.0;
  if (inner) {
    MItem* it = knot->list;
    for (i = outer/2; i < inner + outer/2; i++, it = it->next) {
      if (it->n >= 0) error("curve #%d: invalid knot vector item.", id);
      nurbs->kv[i] = it->val;
    }
  }
  for (i = outer/2 + inner; i < nurbs->nk; i++)
    nurbs->kv[i] = 1.0;

  nurbs->ref = 0;
  return nurbs;
}


//// mesh::create //////////////////////////////////////////////////////////////////////////////////

void Mesh::create(int nv, double2* verts, int nt, int4* tris,
                  int nq, int5* quads, int nm, int3* mark)
{
  free();

  // initialize hash table
  int size = 16;
  while (size < 2*nv) size *= 2;
  HashTable::init(size);

  // create vertex nodes
  for (int i = 0; i < nv; i++)
  {
    Node* node = nodes.add();
    assert(node->id == i);
    node->ref = TOP_LEVEL_REF;
    node->type = TYPE_VERTEX;
    node->bnd = 0;
    node->p1 = node->p2 = -1;
    node->next_hash = NULL;
    node->x = verts[i][0];
    node->y = verts[i][1];
  }
  ntopvert = nv;

  // create triangles
  for (int i = 0; i < nt; i++)
    create_triangle(tris[i][3], &nodes[tris[i][0]], &nodes[tris[i][1]], &nodes[tris[i][2]], NULL);

  // create quads
  for (int i = 0; i < nq; i++)
    create_quad(quads[i][4], &nodes[quads[i][0]], &nodes[quads[i][1]], &nodes[quads[i][2]], &nodes[quads[i][3]], NULL);

  // set boundary markers
  for (int i = 0; i < nm; i++)
  {
    Node* en = peek_edge_node(mark[i][0], mark[i][1]);
    if (en == NULL) error("boundary data error (edge does not exist)");
    en->marker = mark[i][2];

    if (en->marker > 0)
    {
      nodes[mark[i][0]].bnd = 1;
      nodes[mark[i][1]].bnd = 1;
      en->bnd = 1;
    }
  }

  nbase = nactive = ninitial = nt + nq;
  seq = g_mesh_seq++;
}


//// mesh::save ////////////////////////////////////////////////////////////////////////////////////

void Mesh::save_refinements(FILE* f, Element* e, int id, bool& first)
{
  if (e->active) return;
  fprintf(f, first ? "refinements =\n{\n" : ",\n"); first = false;
  if (e->bsplit())
  {
    fprintf(f, "  { %d, 0 }", id);
    int sid = seq; seq += 4;
    for (int i = 0; i < 4; i++)
      save_refinements(f, e->sons[i], sid+i, first);
  }
  else if (e->hsplit())
  {
    fprintf(f, "  { %d, 1 }", id);
    int sid = seq; seq += 2;
    save_refinements(f, e->sons[0], sid, first);
    save_refinements(f, e->sons[1], sid+1, first);
  }
  else
  {
    fprintf(f, "  { %d, 2 }", id);
    int sid = seq; seq += 2;
    save_refinements(f, e->sons[2], sid, first);
    save_refinements(f, e->sons[3], sid+1, first);
  }
}


void Mesh::save_nurbs(FILE* f, int p1, int p2, Nurbs* nurbs)
{
  if (nurbs->arc)
  {
    fprintf(f, "  { %d, %d, %.16g }", p1, p2, nurbs->angle);
  }
  else
  {
    int inner = nurbs->np - 2;
    int outer = nurbs->nk - inner;
    fprintf(f, "  { %d, %d, %d, { ", p1, p2, nurbs->degree);
    for (int i = 1; i < nurbs->np-1; i++)
      fprintf(f, "{ %.16g, %.16g, %.16g }%s ",
                 nurbs->pt[i][0], nurbs->pt[i][1], nurbs->pt[i][2],
                 i < nurbs->np-2 ? "," : "");

    fprintf(f, "}, { ", nurbs->nk - 2*(nurbs->degree+1));
    int max = nurbs->nk - (nurbs->degree+1);
    for (int i = nurbs->degree+1; i < max; i++)
      fprintf(f, "%.16g%s", nurbs->kv[i], i < max-1 ? "," : "");
    fprintf(f, "} }");
  }
}


static bool is_twin_nurbs(Element* e, int i)
{
  // on internal edges, where there are two Nurbs', we only save one of them
  return e->cm->nurbs[i]->twin && e->en[i]->ref == 2;
}


void Mesh::save(const char* filename)
{
  int i, mrk;
  Element* e;

  // open output file
  FILE* f = fopen(filename, "w");
  if (f == NULL) error("Could not create mesh file.");
  //fprintf(f, "# hermes2d saved mesh\n\n");

  // save vertices
  fprintf(f, "vertices =\n{\n");
  for (i = 0; i < ntopvert; i++)
    fprintf(f, "  { %.16g, %.16g }%s\n", nodes[i].x, nodes[i].y, (i < ntopvert-1 ? "," : ""));

  // save elements
  fprintf(f, "}\n\nelements =\n{");
  bool first = true;
  for (i = 0; i < get_num_base_elements(); i++)
  {
    const char* nl = first ? "\n" : ",\n";  first = false;
    e = get_element_fast(i);
    if (!e->used)
      fprintf(f, "%s  { }", nl);
    else if (e->is_triangle())
      fprintf(f, "%s  { %d, %d, %d, %d }", nl, e->vn[0]->id, e->vn[1]->id, e->vn[2]->id, e->marker);
    else
      fprintf(f, "%s  { %d, %d, %d, %d, %d }", nl, e->vn[0]->id, e->vn[1]->id, e->vn[2]->id, e->vn[3]->id, e->marker);
  }

  // save boundary markers
  fprintf(f, "\n}\n\nboundaries =\n{");
  first = true;
  for_all_base_elements(e, this)
    for (i = 0; i < e->nvert; i++)
      if ((mrk = get_base_edge_node(e, i)->marker)) {
        const char* nl = first ? "\n" : ",\n";  first = false;
        fprintf(f, "%s  { %d, %d, %d }", nl, e->vn[i]->id, e->vn[e->next_vert(i)]->id, mrk);
      }
  fprintf(f, "\n}\n\n");

  // save curved edges
  first = true;
  for_all_base_elements(e, this)
    if (e->is_curved())
      for (i = 0; i < e->nvert; i++)
        if (e->cm->nurbs[i] != NULL && !is_twin_nurbs(e, i)) {
          fprintf(f, first ? "curves =\n{\n" : ",\n");  first = false;
          save_nurbs(f, e->vn[i]->id, e->vn[e->next_vert(i)]->id, e->cm->nurbs[i]);
        }
  if (!first) fprintf(f, "\n}\n\n");

  // save refinements
  unsigned temp = seq;
  seq = nbase;
  first = true;
  for_all_base_elements(e, this)
    save_refinements(f, e, e->id, first);
  if (!first) fprintf(f, "\n}\n\n");

  seq = temp;
  fclose(f);
}


//// mesh copy /////////////////////////////////////////////////////////////////////////////////////

void Mesh::copy(const Mesh* mesh)
{
  int i;

  free();

  // copy nodes and elements
  HashTable::copy(mesh);
  elements.copy(mesh->elements);

  Element* e;
  for_all_elements(e, this)
  {
    // update vertex node pointers
    for (i = 0; i < e->nvert; i++)
      e->vn[i] = &nodes[e->vn[i]->id];

    if (e->active)
    {
      // update edge node pointers
      for (i = 0; i < e->nvert; i++)
        e->en[i] = &nodes[e->en[i]->id];
    }
    else
    {
      // update son pointers
      for (i = 0; i < 4; i++)
        if (e->sons[i] != NULL)
          e->sons[i] = &elements[e->sons[i]->id];
    }

    // copy CurvMap, update its parent
    if (e->cm != NULL)
    {
      e->cm = new CurvMap(e->cm);
      if (!e->cm->toplevel)
        e->cm->parent = &elements[e->cm->parent->id];
    }
  }

  // update element pointers in edge nodes
  Node* node;
  for_all_edge_nodes(node, this)
    for (i = 0; i < 2; i++)
      if (node->elem[i] != NULL)
        node->elem[i] = &elements[node->elem[i]->id];

  nbase = mesh->nbase;
  nactive = mesh->nactive;
  ntopvert = mesh->ntopvert;
  ninitial = mesh->ninitial;
  seq = mesh->seq;
}


Node* Mesh::get_base_edge_node(Element* base, int edge)
{
  while (!base->active) // we need to go down to an active element
  {
    int son1, son2;
    get_edge_sons(base, edge, son1, son2);
    base = base->sons[son1];
  }
  return base->en[edge];
}


void Mesh::copy_base(Mesh* mesh)
{
  free();
  HashTable::init();

  // copy top-level vertex nodes
  for (int i = 0; i < mesh->get_max_node_id(); i++)
  {
    Node* node = &(mesh->nodes[i]);
    if (node->ref < TOP_LEVEL_REF) break;
    Node* newnode = nodes.add();
    assert(newnode->id == i && node->type == TYPE_VERTEX);
    memcpy(newnode, node, sizeof(Node));
    newnode->ref = TOP_LEVEL_REF;
  }

  // copy base elements
  Element* e;
  for_all_base_elements(e, mesh)
  {
    Element* enew;
    Node *v0 = &nodes[e->vn[0]->id], *v1 = &nodes[e->vn[1]->id], *v2 = &nodes[e->vn[2]->id];
    if (e->is_triangle())
      enew = create_triangle(e->marker, v0, v1, v2, NULL);
    else
      enew = create_quad(e->marker, v0, v1, v2, &nodes[e->vn[3]->id], NULL);

    // copy edge markers
    for (int j = 0; j < e->nvert; j++)
    {
      Node* en = get_base_edge_node(e, j);
      enew->en[j]->bnd = en->bnd; // copy bnd data from the active el.
      enew->en[j]->marker = en->marker;
    }

    enew->userdata = e->userdata;
    if (e->is_curved())
      enew->cm = new CurvMap(e->cm);
  }

  nbase = nactive = ninitial = mesh->nbase;
  ntopvert = mesh->ntopvert;
  seq = g_mesh_seq++;
}


//// save_raw, load_raw ////////////////////////////////////////////////////////////////////////////

void Mesh::save_raw(FILE* f)
{
  int i, nn, mm;
  int null = -1;

  assert(sizeof(int) == 4);
  assert(sizeof(double) == 8);

  hermes2d_fwrite("H2DM\001\000\000\000", 1, 8, f);

  #define output(n, type) \
    hermes2d_fwrite(&(n), sizeof(type), 1, f)

  output(nbase, int);
  output(ntopvert, int);
  output(nactive, int);

  nn = nodes.get_num_items();
  mm = nodes.get_size();
  output(nn, int);
  output(mm, int);

  // dump all nodes
  Node* n;
  for_all_nodes(n, this)
  {
    output(n->id, int);
    unsigned bits = n->ref | (n->type << 29) | (n->bnd << 30) | (n->used << 31);
    output(bits, unsigned);

    if (n->type == TYPE_VERTEX)
    {
      output(n->x, double);
      output(n->y, double);
    }
    else
    {
      output(n->marker, int);
      output(n->elem[0] ? n->elem[0]->id : null, int);
      output(n->elem[1] ? n->elem[1]->id : null, int);
    }

    output(n->p1, int);
    output(n->p2, int);
  }

  nn = elements.get_num_items();
  mm = elements.get_size();
  output(nn, int);
  output(mm, int);

  // dump all elements
  Element* e;
  for (int id = 0; id < get_max_element_id(); id++)
  {
    if ((e = get_element_fast(id))->used || id < nbase)
    {
      output(e->id, int);
      unsigned bits = e->nvert | (e->active << 30) | (e->used << 31);
      output(bits, unsigned);

      if (e->used)
      {
        output(e->marker, int);
        output(e->userdata, int);
        output(e->iro_cache, int);

        for (i = 0; i < e->nvert; i++)
          output(e->vn[i]->id, int);

        if (e->active)
          for (i = 0; i < e->nvert; i++)
            output(e->en[i]->id, int);
        else
          for (i = 0; i < 4; i++)
            output(e->sons[i] ? e->sons[i]->id : null, int);

        if (e->is_curved()) error("Not implemented for curved elements yet.");
      }
    }
  }

  // TODO: curved elements

  #undef output
}


void Mesh::load_raw(FILE* f)
{
  int i, j, nv, mv, ne, me, id;

  assert(sizeof(int) == 4);
  assert(sizeof(double) == 8);

  // check header
  struct { char magic[4]; int ver; } hdr;
  hermes2d_fread(&hdr, sizeof(hdr), 1, f);
  if (hdr.magic[0] != 'H' || hdr.magic[1] != '2' || hdr.magic[2] != 'D' || hdr.magic[3] != 'M')
    error("Not a Hermes2D raw mesh file.");
  if (hdr.ver > 1)
    error("Unsupported file version.");

  #define input(n, type) \
    hermes2d_fread(&(n), sizeof(type), 1, f)

  free();

  input(nbase, int);
  input(ntopvert, int);
  input(nactive, int);

  input(nv, int);
  input(mv, int);
  nodes.force_size(mv);

  // load nodes
  for (i = 0; i < nv; i++)
  {
    input(id, int);
    if (id < 0 || id >= mv) error("Corrupt data.");
    Node* n = &(nodes[id]);
    n->id = id;
    n->used = 1;

    unsigned bits;
    input(bits, unsigned);
    n->ref  =  bits & 0x1fffffff;
    n->type = (bits >> 29) & 0x1;
    n->bnd  = (bits >> 30) & 0x1;

    if (n->type == TYPE_VERTEX)
    {
      input(n->x, double);
      input(n->y, double);
    }
    else
    {
      input(n->marker, int);
      n->elem[0] = n->elem[1] = NULL;
      input(n->elem[0], int);
      input(n->elem[1], int);
    }

    input(n->p1, int);
    input(n->p2, int);
  }
  nodes.post_load_scan();

  int hsize = DEFAULT_HASH_SIZE;
  while (hsize < nv) hsize *= 2;
  HashTable::init(hsize);
  HashTable::rebuild();

  input(ne, int);
  input(me, int);
  elements.force_size(me);

  // load elements
  for (i = 0; i < ne; i++)
  {
    input(id, int);
    if (id < 0 || id >= me) error("Corrupt data.");
    Element* e = &(elements[id]);
    e->id = id;

    unsigned bits;
    input(bits, unsigned);
    e->nvert  =  bits & 0x3fffffff;
    e->active = (bits >> 30) & 0x1;
    e->used   = (bits >> 31) & 0x1;

    if (e->used)
    {
      input(e->marker, int);
      input(e->userdata, int);
      input(e->iro_cache, int);

      // load vertex node ids
      for (j = 0; j < e->nvert; j++)
      {
        input(id, int);
        if (id < 0 || id >= mv) error("Corrupt data.");
        e->vn[j] = get_node(id);
      }

      if (e->active)
      {
        // load edge node ids
        for (j = 0; j < e->nvert; j++)
        {
          input(id, int);
          if (id < 0 || id >= mv) error("Corrupt data.");
          e->en[j] = get_node(id);
        }
      }
      else
      {
        // load son ids
        for (j = 0; j < 4; j++)
        {
          input(id, int);
          if (id < 0)
            e->sons[j] = NULL;
          else if (id < me)
            e->sons[j] = &(elements[id]);
          else
            error("Corrupt data.");
        }
      }
    }
  }
  elements.post_load_scan(nbase);

  // update edge node element pointers
  Node* n;
  for_all_edge_nodes(n, this)
    for (j = 0; j < 2; j++)
      if ((int) (long) n->elem[j] == -1)
        n->elem[j] = NULL;
      else
        n->elem[j] = get_element((int) (long) n->elem[j]);

  #undef input
  seq++;
}


//// free //////////////////////////////////////////////////////////////////////////////////////////

void Mesh::free()
{
  Element* e;
  for_all_elements(e, this)
    if (e->cm != NULL)
    {
      delete e->cm;
      e->cm = NULL; // fixme!!!
    }

  elements.free();
  HashTable::free();
}

void Mesh::copy_refine(Mesh* mesh)
{
  free();
  HashTable::copy(mesh);

  // copy base elements
  Element* e;

  for_all_refine_elements(e, mesh)
  {
    Element* enew;
    Node *v0 = &nodes[e->vn[0]->id], *v1 = &nodes[e->vn[1]->id], *v2 = &nodes[e->vn[2]->id];
    if (e->is_triangle())
      enew = create_triangle(e->marker, v0, v1, v2, NULL);
    else
      enew = create_quad(e->marker, v0, v1, v2, &nodes[e->vn[3]->id], NULL);

    // copy edge markers
    for (int j = 0; j < e->nvert; j++)
    {
      Node* en = get_base_edge_node(e, j);
      enew->en[j]->bnd = en->bnd;
      enew->en[j]->marker = en->marker;
    }

    enew->userdata = e->userdata;
    if (e->is_curved())
      enew->cm = new CurvMap(e->cm);
  }

  nbase = nactive = ninitial = mesh->nbase;
  ntopvert = mesh->ntopvert;
  seq = g_mesh_seq++;
}

////convert a triangle element into three quadrilateral elements///////

void Mesh::convert_to_quads(int refinement)
{
  int i;
  Element* e;
  int temp_n = 900000;
  int temp_compare[temp_n];
  int temp_exist[temp_n];
  int temp_i = 0;
  int temp_j = 0;
  int temp_k = 0;
  int temp_count = 0;
  int max_node_num = 0;

  Mesh tmp;
  memset(temp_compare, 0, sizeof(temp_compare));
  memset(temp_exist, 0, sizeof(temp_exist));

  for (temp_i = 0; temp_i < temp_n; temp_i++)
  {
    temp_compare[temp_i] = temp_i;
  }

  elements.set_append_only(true);
  for_all_active_elements(e, this)
    refine_element_to_quads(e->id, refinement);
  elements.set_append_only(false);

  tmp.copy_refine(this);
  copy(&tmp);
  // get all nodes number and value
  for (i = get_num_base_elements(); i < get_max_element_id(); i++)
  {
    e = get_element_fast(i);
    if (e->is_triangle())
           temp_k = 3;
    else
           temp_k = 4;

    for (temp_j = 0; temp_j < temp_k; temp_j++)
    {
      if (temp_compare[e->vn[temp_j]->id] == e->vn[temp_j]->id)
        temp_exist[e->vn[temp_j]->id] = e->vn[temp_j]->id;
    }
  }

  for (i = 1; i < temp_n; i++)
  {
    if (temp_exist[i] != 0)
      temp_count++;
  }

  // get the total count of node number
  for (i = 1; i < temp_n; i++)
  {
    if (temp_exist[i] > temp_exist[i-1])
    {
        max_node_num = temp_exist[i] + 1;
    }
  }
  temp_count += 1;
  ntopvert = max_node_num;
  nbase = get_max_element_id();
}

////convert a quad element into two triangle elements///////
void Mesh::convert_to_triangles()
{
  int i;
  Element* e;
  int temp_n = 900000;
  int temp_compare[temp_n];
  int temp_exist[temp_n];
  int temp_i = 0;
  int temp_j = 0;
  int temp_k = 0;
  int temp_count = 0;
  int max_node_num = 0;

  Mesh tmp;
  memset(temp_compare, 0, sizeof(temp_compare));
  memset(temp_exist, 0, sizeof(temp_exist));

  for (temp_i = 0; temp_i < temp_n; temp_i++)
  {
    temp_compare[temp_i] = temp_i;
  }

  elements.set_append_only(true);
  for_all_active_elements(e, this)
    refine_element_to_triangles(e->id);
  elements.set_append_only(false);

  tmp.copy_refine(this);
  copy(&tmp);
  // get all nodes number and value
  for (i = get_num_base_elements(); i < get_max_element_id(); i++)
  {
    e = get_element_fast(i);
    if (e->is_triangle())
           temp_k = 3;
    else
           temp_k = 4;

    for (temp_j = 0; temp_j < temp_k; temp_j++)
    {
      if (temp_compare[e->vn[temp_j]->id] == e->vn[temp_j]->id)
        temp_exist[e->vn[temp_j]->id] = e->vn[temp_j]->id;
    }
  }

  for (i = 1; i < temp_n; i++)
  {
    if (temp_exist[i] != 0)
      temp_count++;
  }

  // get the total count of node number
  for (i = 1; i < temp_n; i++)
  {
    if (temp_exist[i] > temp_exist[i-1])
    {
        max_node_num = temp_exist[i] + 1;
    }
  }
  temp_count += 1;
  ntopvert = max_node_num;
  nbase = get_max_element_id();
}

