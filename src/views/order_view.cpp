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
#include "../common.h"
#include "order_view.h"

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
  center_mesh(ord.get_vertices(), ord.get_num_vertices());
  ord.unlock_data();

  create();
  init_order_palette();
  update_layout();
  refresh();
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
    glColor3f(0.4f, 0.4f, 0.4f);
  else if (pal_type == 1)
    glColor3f(1.0f, 1.0f, 1.0f);
  else
    glColor3f(0.0f, 0.0f, 0.0f);
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
      refresh();
      break;
    
    case 'm':
      b_orders = !b_orders;
      refresh();
      break;
        
    case 'p':
    {
      pal_type++;
      if (pal_type > 2) pal_type = 0; 
      init_order_palette();
      refresh();
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
  ord.lock_data();
  center_mesh(ord.get_vertices(), ord.get_num_vertices());
  ord.unlock_data();
  
  create();
  init_order_palette();
  update_layout();
  refresh();
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

#endif // NOGLUT
