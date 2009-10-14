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

// $Id: view4.cpp 1086 2008-10-21 09:05:44Z jakub $

#ifndef NOGLUT

#include <GL/freeglut.h>
#include "common.h"
#include "view.h"
#include "discrete.h"


//// MeshView //////////////////////////////////////////////////////////////////////////////////////

MeshView::MeshView(const char* title, int x, int y, int width, int height)
        : View(title, x, y, width, height)
{
  nodes = elems = NULL;
  b_scale = false;
  b_ids = false;
  b_markers = true;
}


void MeshView::show(Mesh* mesh)
{
  glut_init();
  update_layout();

  Solution sln;
  sln.set_zero(mesh);
  lin.process_solution(&sln);
  lin.lock_data();
  center_mesh(lin.get_vertices(), lin.get_num_vertices());
  lin.unlock_data();

  int i;

  if (elems != NULL) delete [] elems;
  ne = mesh->get_max_element_id()+1;
  elems = new ObjInfo[ne];
  for (i = 0; i < ne; i++)
    elems[i].id = -1;

  Element* e;
  for_all_active_elements(e, mesh)
  {
    ObjInfo* oi = elems + e->id;
    oi->id = e->id;
    oi->type = e->marker;
    oi->x = oi->y = 0.0;
    for (i = 0; i < e->nvert; i++) {
      oi->x += e->vn[i]->x;
      oi->y += e->vn[i]->y;
    }
    oi->x /= e->nvert;
    oi->y /= e->nvert;
  }

  create();
  wait_for_draw();
}


void MeshView::on_display()
{
  set_ortho_projection();
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  // transform all vertices
  lin.lock_data();
  int i, nv = lin.get_num_vertices();
  double3* vert = lin.get_vertices();
  double2* tvert = new double2[nv];
  for (i = 0; i < nv; i++)
  {
    tvert[i][0] = transform_x(vert[i][0]);
    tvert[i][1] = transform_y(vert[i][1]);
  }

  // draw all triangles
  int3* tris = lin.get_triangles();
  glColor3f(0.9, 0.9, 0.9);
  glBegin(GL_TRIANGLES);
  for (i = 0; i < lin.get_num_triangles(); i++)
  {
    glVertex2d(tvert[tris[i][0]][0], tvert[tris[i][0]][1]);
    glVertex2d(tvert[tris[i][1]][0], tvert[tris[i][1]][1]);
    glVertex2d(tvert[tris[i][2]][0], tvert[tris[i][2]][1]);
  }
  glEnd();

  // draw all edges
  glLineStipple(5, 0x5555);
  int3* edges = lin.get_edges();
  for (i = 0; i < lin.get_num_edges(); i++)
  {
    int mrk = b_markers ? edges[i][2] : 0;

    if (!edges[i][2] &&
        (tvert[edges[i][0]][1] == tvert[edges[i][1]][1] && tvert[edges[i][0]][0] < tvert[edges[i][1]][0] ||
         tvert[edges[i][0]][1] < tvert[edges[i][1]][1])) continue;

    float* color = get_marker_color(mrk);
    glColor3f(color[0], color[1], color[2]);
    glLineWidth(mrk ? 1.5 : 1.0);
    glBegin(GL_LINES);
      glVertex2d(tvert[edges[i][0]][0], tvert[edges[i][0]][1]);
      glVertex2d(tvert[edges[i][1]][0], tvert[edges[i][1]][1]);
    glEnd();

    if (mrk)
    {
      glEnable(GL_LINE_STIPPLE);
      glColor3f(0.4, 0.4, 0.4);
      glBegin(GL_LINES);
        glVertex2d(tvert[edges[i][0]][0], tvert[edges[i][0]][1]);
        glVertex2d(tvert[edges[i][1]][0], tvert[edges[i][1]][1]);
      glEnd();
      glDisable(GL_LINE_STIPPLE);
    }
  }
  glLineWidth(1.0);

  // draw element ids
  if (b_ids)
  {
    glColor3f(0, 0, 0);
    for (i = 0; i < ne; i++)
    {
      if (elems[i].id < 0) continue;
      char text[20];
      sprintf(text, "#%d", elems[i].id);
      draw_text(transform_x(elems[i].x), transform_y(elems[i].y), text, 0);
    }
  }

  delete [] tvert;
  lin.unlock_data();
}


void MeshView::on_key_down(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 'c':
      lin.lock_data();
      center_mesh(lin.get_vertices(), lin.get_num_vertices());
      lin.unlock_data();
      post_redisplay();
      //reset_zoom();
      break;

    case 'b':
      b_markers = !b_markers;
      post_redisplay();
      break;

    case 'm':
      b_ids = !b_ids;
      post_redisplay();
      break;

    default:
      View::on_key_down(key, x, y);
      break;
  }
}


float* MeshView::get_marker_color(int marker)
{
  static float edgecol[3] = { 0.3, 0.3, 0.3 };
  static float randcol[3];
  static float mc[8][3] =
  {
    { 1.0, 0.3, 0.3 },
    { 0.0, 0.9, 0.0 },
    { 0.0, 0.0, 0.7 },
    { 1.0, 1.0, 0.2 },
    { 0.7, 0.0, 0.0 },
    { 0.0, 0.5, 0.0 },
    { 0.3, 0.5, 1.0 },
    { 0.8, 0.8, 0.0 },
  };

  if (marker == 0)
    return edgecol;
  else if (marker > 0 && marker < 8)
    return mc[marker];
  else {
    srand(marker+2);
    randcol[0] = (float) rand() / RAND_MAX;
    randcol[1] = (float) rand() / RAND_MAX;
    randcol[2] = (float) rand() / RAND_MAX;
    return randcol;
  }
}


MeshView::~MeshView()
{
  if (nodes != NULL) delete [] nodes;
  if (elems != NULL) delete [] elems;
}


const char* MeshView::get_help_text() const
{
  return
  "MeshView\n"
  "Controls:\n"
  "  Left mouse - pan\n"
  "  Right mouse - zoom\n"
  "  C - center image\n"
  "  H - render high-quality frame\n"
  "  B - toggle boundary markers\n"
  "  M - toggle element IDs\n"
  "  S - save screenshot\n"
  "  F1 - this help\n"
  "  Esc, Q - quit";
}


//// BaseView //////////////////////////////////////////////////////////////////////////////////////

BaseView::BaseView(const char* title, int x, int y, int width, int height)
        : ScalarView(title, x, y, width, height)
{
  pss = NULL;
  sln = NULL;
  lines = true;
}


void BaseView::show(Space* space, double eps, int item)
{
  free();
  // todo: copy space
  pss = new PrecalcShapeset(space->get_shapeset());
  sln = new Solution();
  ndofs = space->get_num_dofs();
  base_index = 0;
  this->space = space;
  this->eps = eps;
  this->item = item;
  update_solution();
}


void BaseView::free()
{
  if (pss != NULL) { delete pss; pss = NULL; }
  if (sln != NULL) { delete sln; sln = NULL; }
}


void BaseView::update_solution()
{
  scalar* vec = new scalar[ndofs];
  memset(vec, 0, sizeof(scalar) * ndofs);
  if (base_index >= 0)
  {
    if (base_index < ndofs) vec[base_index] = 1.0;
    sln->set_fe_solution(space, pss, vec, 0.0);
  }
  else
  {
    sln->set_fe_solution(space, pss, vec, 1.0);
  }

  ScalarView::show(sln, eps, item);
  update_title();
}


void BaseView::update_title()
{
  char text[500];
  sprintf(text, "%s - dof = %d%s", title.c_str(), base_index,
          (base_index < 0) ? " (Dirichlet lift)" : "");
  set_title_internal(text);
}


void BaseView::on_special_key(int key, int x, int y)
{
  switch (key)
  {
    case GLUT_KEY_LEFT:
      if (base_index > -1) base_index--;
      update_solution();
      break;

    case GLUT_KEY_RIGHT:
      if (base_index < ndofs-1) base_index++;
      update_solution();
      break;

    default:
      ScalarView::on_special_key(key, x, y);
  }
}


const char* BaseView::get_help_text() const
{
  return
  "BaseView\n\n"
  "Controls:\n"
  "  Left mouse - pan\n"
  "  Right mouse - zoom\n"
  "  Left arrow - previous basis function\n"
  "  Right arrow - next basis function\n"
  "  3 - toggle 3D mode\n"
  "  C - center image\n"
  "  F - toggle smooth palette\n"
  "  H - render high-quality frame\n"
  "  M - toggle mesh\n"
  "  P - cycle palettes\n"
  "  S - save screenshot\n"
  "  F1 - this help\n"
  "  Esc, Q - quit\n\n"
  "3D mode:\n"
  "  Left mouse - rotate\n"
  "  Right mouse - zoom\n"
  "  Middle mouse - pan\n"
  "  * - increase Z scale\n"
  "  / - decrease Z scale";
}


//// VecBaseView //////////////////////////////////////////////////////////////////////////////////////

void VectorBaseView::show(Space* space)
{
  free();
  pss = new PrecalcShapeset(space->get_shapeset());
  sln = new Solution();
  this->space = space;
  ndofs = space->get_num_dofs();
  base_index = 0;
  update_solution();
}


void VectorBaseView::free()
{
  if (pss != NULL) { delete [] pss; pss = NULL; }
  if (sln != NULL) { delete [] sln; sln = NULL; }
}


void VectorBaseView::update_solution()
{
  scalar* vec = new scalar[ndofs + 1];
  memset(vec, 0, sizeof(scalar) * (ndofs + 1));
  if (base_index >= -1 && base_index < ndofs)
    vec[base_index + 1] = 1.0;
  sln->set_fe_solution(space, pss, vec);

  VectorView::show(sln,  sln, 0.001, FN_VAL_0, FN_VAL_1);
  update_title();
}


void VectorBaseView::update_title()
{
  char text[500];
  sprintf(text, "%s - dof = %d%s", title.c_str(), base_index,
          (base_index < 0) ? " (Dirichlet lift)" : "");
  set_title_internal(text);
}


void VectorBaseView::on_special_key(int key, int x, int y)
{
  switch (key)
  {
    case GLUT_KEY_LEFT:
      if (base_index > -1) base_index--;
      update_solution();
      break;

    case GLUT_KEY_RIGHT:
      if (base_index < ndofs-1) base_index++;
      update_solution();
      break;

    default:
      VectorView::on_special_key(key, x, y);
  }
}


const char* VectorBaseView::get_help_text() const
{
  return
  "VectorBaseView\n\n"
  "Controls:\n"
  "  Left mouse - pan\n"
  "  Right mouse - zoom\n"
  "  Left arrow - previous basis function\n"
  "  Right arrow - next basis function\n"
  "  C - center image\n"
  "  F - toggle smooth palette\n"
  "  X - toggle hexagonal grid\n"
  "  H - render high-quality frame\n"
  "  M - toggle mesh\n"
  "  P - cycle palettes\n"
  "  S - save screenshot\n"
  "  F1 - this help\n"
  "  Esc, Q - quit";
}


//// OrderView /////////////////////////////////////////////////////////////////////////////////////

OrderView::OrderView(const char* title, int x, int y, int width, int height)
         : View(title, x, y, width, height)
{
  b_scale = true;
  b_orders = false;
  scale_width = 36;
  scale_box_height = 25;
  scale_box_skip = 9;
}


static int order_palette[] =
{
  0x000000, 0x000684, 0x3250fc, 0x36c4ee, 0x04eabc,
  0x62ff2a, 0xfdff07, 0xffa044, 0xff1111, 0xb02c2c, 0x820f97
};


void OrderView::show(Space* space)
{
  if (!space->is_up_to_date())
    error("The space is not up to date.");

  ord.lock_data();
  ord.process_solution(space);
  init_order_palette();
  glut_init();
  update_layout();

  if (window_id < 0)
  {
    center_mesh(ord.get_vertices(), ord.get_num_vertices());
  }

  create();
  ord.unlock_data();
  wait_for_draw();
}


void OrderView::init_order_palette()
{
  int min = 1, max = 1;
  double3* vert = ord.get_vertices();
  for (int i = 0; i < ord.get_num_vertices(); i++)
  {
    if ((int) vert[i][2] < min) min = (int) vert[i][2];
    if ((int) vert[i][2] > max) max = (int) vert[i][2];
  }

  num_boxes = max - min + 1;
  char* buf = text_buffer;
  for (int i = 0; i < num_boxes; i++)
  {
    //if (pal_type != 0)
      memcpy(order_colors[i+min], get_palette_color((i + min - 1) / 9.0), 3*sizeof(float));
    /*else {
      order_colors[i+min][0] = (float) (order_palette[i+min] >> 16) / 0xff;
      order_colors[i+min][1] = (float) ((order_palette[i+min] >> 8) & 0xff) / 0xff;
      order_colors[i+min][2] = (float) (order_palette[i+min] & 0xff) / 0xff;
    }*/

    sprintf(buf, "%d", i + min);
    box_names[i] = buf;
    buf += strlen(buf)+1;
  }

  scale_height = num_boxes * scale_box_height + (num_boxes-1) * scale_box_skip;
  order_min = min;
}


void OrderView::on_display()
{
  set_ortho_projection();
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  // transform all vertices
  ord.lock_data();
  int i, nv = ord.get_num_vertices();
  double3* vert = ord.get_vertices();
  double2* tvert = new double2[nv];
  for (i = 0; i < nv; i++)
  {
    tvert[i][0] = transform_x(vert[i][0]);
    tvert[i][1] = transform_y(vert[i][1]);
  }

  // draw all triangles
  int3* tris = ord.get_triangles();
  glBegin(GL_TRIANGLES);
  const float* color;
  for (i = 0; i < ord.get_num_triangles(); i++)
  {
    color = order_colors[(int) vert[tris[i][0]][2]];
    //get_palette_color((vert[tris[i][0]][2] - 1) / 9.0);
    glColor3f(color[0], color[1], color[2]);

    glVertex2d(tvert[tris[i][0]][0], tvert[tris[i][0]][1]);
    glVertex2d(tvert[tris[i][1]][0], tvert[tris[i][1]][1]);
    glVertex2d(tvert[tris[i][2]][0], tvert[tris[i][2]][1]);
  }
  glEnd();

  // draw all edges
  if (pal_type == 0)
    glColor3f(0.4, 0.4, 0.4);
  else if (pal_type == 1)
    glColor3f(1.0, 1.0, 1.0);
  else
    glColor3f(0.0, 0.0, 0.0);
  glBegin(GL_LINES);
  int3* edges = ord.get_edges();
  for (i = 0; i < ord.get_num_edges(); i++)
  {
    glVertex2d(tvert[edges[i][0]][0], tvert[edges[i][0]][1]);
    glVertex2d(tvert[edges[i][1]][0], tvert[edges[i][1]][1]);
  }
  glEnd();

  // draw labels
  if (b_orders)
  {
    int* lvert;
    char** ltext;
    double2* lbox;
    int nl = ord.get_labels(lvert, ltext, lbox);
    for (i = 0; i < nl; i++)
      if (lbox[i][0] * scale > get_text_width(ltext[i]) &&
          lbox[i][1] * scale > 13)
      {
        //color = get_palette_color((vert[lvert[i]][2] - 1) / 9.0);
        color = order_colors[(int) vert[lvert[i]][2]];
        if ((color[0] + color[1] + color[2]) / 3 > 0.5)
          glColor3f(0, 0, 0);
        else
          glColor3f(1, 1, 1);

        draw_text(tvert[lvert[i]][0], tvert[lvert[i]][1], ltext[i], 0);
      }
  }

  delete [] tvert;
  ord.unlock_data();
}


int OrderView::measure_scale_labels()
{
  return 0;
}

void OrderView::scale_dispatch()
{
  draw_discrete_scale(num_boxes, box_names, order_colors + order_min);
}


void OrderView::on_key_down(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 'c':
      ord.lock_data();
      center_mesh(ord.get_vertices(), ord.get_num_vertices());
      ord.unlock_data();
      post_redisplay();
      break;

    case 'm':
      b_orders = !b_orders;
      post_redisplay();
      break;

    case 'p':
    {
      pal_type++;
      if (pal_type > 2) pal_type = 0;
      init_order_palette();
      post_redisplay();
      break;
    }

    default:
      View::on_key_down(key, x, y);
      break;
  }
}


void OrderView::load_data(const char* filename)
{
  ord.load_data(filename);
  init_order_palette();
  glut_init();
  update_layout();

  if (window_id < 0)
  {
    ord.lock_data();
    center_mesh(ord.get_vertices(), ord.get_num_vertices());
    ord.unlock_data();
  }

  create();
  wait_for_draw();
}


void OrderView::save_data(const char* filename)
{
  if (ord.get_num_triangles() <= 0)
    error("No data to save.");
  ord.save_data(filename);
}


void OrderView::save_numbered(const char* format, int number)
{
  char buffer[1000];
  sprintf(buffer, format, number);
  save_data(buffer);
}


const char* OrderView::get_help_text() const
{
  return
  "OrderView\n\n"
  "Controls:\n"
  "  Left mouse - pan\n"
  "  Right mouse - zoom\n"
  "  C - center image\n"
  "  M - toggle element orders\n"
  "  H - render high-quality frame\n"
  "  P - cycle palettes\n"
  "  S - save screenshot\n"
  "  F1 - this help\n"
  "  Esc, Q - quit";
}


//// Matrix View ///////////////////////////////////////////////////////////////////////////////////

/*MatrixView::MatrixView(const char* title, int x, int y, int width, int height)
          : View(title, x, y, width, height)
{
  Ap = NULL;
  Ai = NULL;
  Ax = NULL;
  size = 0;
  eps = 1e-12;
  max = 0.0;
  pointsize = 1.0;
  b_scale = false;
}


static inline bool is_zero(scalar x, double eps)
{
#ifndef COMPLEX
  return (fabs(x) < eps);
#else
  return (norm(x) < sqr(eps));
#endif
}

const unsigned MARK_ZERO = (0 << 30),
               MARK_SYMM = (1 << 30),
               MARK_ANTISYMM = ((unsigned) 2 << 30),
               MARK_NOSYMM = ((unsigned) 3 << 30),
               MARK_MASK = (0x3FFFFFFF);


static inline int find_entry(unsigned *Ai, int i, int lo, int hi)
{
  int k;
  while(1)
  {
    k = (lo + hi) >> 1;
    if (i < (Ai[k] & MARK_MASK))
      hi = k - 1;
    else if (i > (Ai[k] & MARK_MASK))
      lo = k + 1;
    else
      return k;
    if (lo > hi)
      return -1;
  }
}


void MatrixView::find_symmetry(int nz, int *sz)
{
  *sz = 0;

  int i,j,k;
  for (i = 0; i < size; i++)
  {
    for (j = 0; j < Ap[i+1] - Ap[i]; j++)
    {
      int row = (Ai[Ap[i] + j] & MARK_MASK);
      if (is_zero(Ax[Ap[i] + j], max*eps))
      {
        Ai[Ap[i] + j] |= MARK_ZERO; // zero entry
        (*sz)++;
      }
      else
      {
        k = find_entry(Ai, i, Ap[row], Ap[row + 1] - 1);
        if (k == -1)
        {
          Ai[Ap[i] + j] |= MARK_NOSYMM; // nonsymmetric entry (other one is not in the matrix at all)
          continue;
        }
        if (is_zero(Ax[k] - Ax[Ap[i] + j], eps))
          Ai[Ap[i] + j] |= MARK_SYMM;  //symmetric entries
        else
        {
          if (is_zero(Ax[k] + Ax[Ap[i] + j],eps))
            Ai[Ap[i] + j] |= MARK_ANTISYMM;   //antisymetric entries
          else
            Ai[Ap[i] + j] |= MARK_NOSYMM;  //nonsymmetric entries (different values)
        }
      }
    }
  }
}


void MatrixView::show(DiscreteProblem *ep)
{
  int *App, *Aii;
  ep->get_matrix(App, Aii, Ax, size);

  int nz = App[size];
  int sz;

  if (Ap != NULL) delete [] Ap;
  Ap = new int[size+1];
  memcpy(Ap, App, (size+1)*sizeof(int));

  if (Ai != NULL) delete [] Ai;
  Ai = new unsigned[nz];
  memcpy(Ai, Aii, nz*sizeof(unsigned));

  max = 0.0;
  for (int i = 0; i < size; i++)
    #ifndef COMPLEX
      if (fabs(Ax[Ap[i]]) > max)
        max = fabs(Ax[Ap[i]]);
    #else
      if (std::abs(Ax[Ap[i]]) > max)
        max = std::abs(Ax[Ap[i]]);
    #endif

  find_symmetry(nz, &sz);

  verts[0][0] = -0.5*pointsize;            verts[0][1] = -0.5*pointsize;
  verts[1][0] = -0.5*pointsize;            verts[1][1] =  0.5*pointsize + size - 1;
  verts[2][0] =  0.5*pointsize + size - 1; verts[2][1] =  0.5*pointsize + size - 1;
  verts[3][0] =  0.5*pointsize + size - 1; verts[3][1] = -0.5*pointsize;

  glut_init();
  update_layout();
  center_mesh(verts, 4);

  create();
  update_title(nz, sz);
  wait_for_draw();
}


void MatrixView::assign_color(int i, int j)
{
  static double c[4][3] = { {0.8, 0.8, 0.8}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 0.0, 0.0} };
  int col = (Ai[Ap[i] + j] >> 30);
  glColor3f(c[col][0], c[col][1], c[col][2]);
}


void MatrixView::on_display()
{
  set_ortho_projection();
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  int i, j;

  for (i = 0; i < size; i++)
  {
    for (j = 0; j < Ap[i+1] - Ap[i]; j++)
    {
      assign_color(i, j);
      int idx = (Ai[Ap[i] + j] & 0x3FFFFFFF);
      double ps = pointsize*scale;
      if (ps < 1.0)
      {
        glPointSize(ps);
        glBegin(GL_POINTS);
        glVertex2d(transform_x(i), transform_y(size -1  - idx));
        glEnd();
      }
      else
      {
        glBegin(GL_QUADS);
        double x = transform_x(i);
        double y = transform_y(size -1  - idx);
        glVertex2d(x - ps/2.0, y - ps/2.0);
        glVertex2d(x + ps/2.0, y - ps/2.0);
        glVertex2d(x + ps/2.0, y + ps/2.0);
        glVertex2d(x - ps/2.0, y + ps/2.0);
        glEnd();
      }
    }
  }

  glBegin(GL_LINE_LOOP);
  glColor3f(0.0, 0.0, 0.0);
  for (i = 0; i < 4; i++)
    glVertex2d(transform_x(verts[i][0]), transform_y(verts[i][1]));
  glEnd();
}


void MatrixView::on_key_down(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 'c':
      center_mesh(verts, 4);
      post_redisplay();
      break;

    default:
      View::on_key_down(key, x, y);
      break;
  }
}


void MatrixView::update_title(int nz, int sz)
{
  char text[500];
  sprintf(text, "%s - matrix size = %d x %d, nonzeros = %d, saved zeros = %d", title.c_str(), size, size, nz, sz);
  set_title_internal(text);
}


MatrixView::~MatrixView()
{
  if (Ai != NULL) delete [] Ai;
  if (Ap != NULL) delete [] Ap;
}


const char* MatrixView::get_help_text() const
{
  return
  "MatrixView\n\n"
  "Controls:\n"
  "  Left mouse - pan\n"
  "  Right mouse - zoom\n"
  "  C - center image\n"
  "  H - render high-quality frame\n"
  "  S - save screenshot\n"
  "  F1 - this help\n"
  "  Esc, Q - quit\n\n"
  "Colors:\n"
  "  Blue - symmetric element\n"
  "  Red - unsymmetric element\n"
  "  Magenta - antisymmetric element\n"
  "  Grey - stored zero element";
}*/


#endif // NOGLUT
