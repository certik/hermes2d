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

#ifndef NOGLUT

#include <GL/freeglut.h>
#include <algorithm>
#include <cmath>
#include <list>
#include "common.h"
#include "view.h"

#ifdef ENABLE_VIEWER_GUI
# include <AntTweakBar.h>
#endif

using namespace std;

//// ScalarView ////////////////////////////////////////////////////////////////////////////////////

ScalarView::ScalarView(const char* title, int x, int y, int width, int height)
           : View(title, x, y, width, height)
           , show_element_info(false), element_id_widget(0)
           , vertex_nodes(0), node_pixel_radius(10), pointed_vertex_node(NULL), pointed_node_widget(0), selected_node_widget(0), node_widget_vert_cnt(32)
#ifdef ENABLE_VIEWER_GUI
           , tw_initialized(false), tw_setup_bar(NULL)
#endif
{
  lines = true;
  pmode = mode3d = false;
  normals = NULL;
  panning = false;

  contours = false;
  cont_orig = 0.0;
  cont_step = 0.2;
  cont_color[0] = 0.0f; cont_color[1] = 0.0f; cont_color[2] = 0.0f;

  show_edges = true;
  edges_color[0] = 0.5f; edges_color[1] = 0.4f; edges_color[2] = 0.4f;

  show_values = true;
}

ScalarView::~ScalarView()
{
  if (normals != NULL) delete [] normals;
  vertex_nodes.clear();

# ifdef ENABLE_VIEWER_GUI
  if (tw_initialized)
  {
    TwDeleteBar(tw_setup_bar);
    TwTerminate();
    tw_initialized = false;
  }
# endif
}

void ScalarView::on_create()
{
# ifdef ENABLE_VIEWER_GUI
  if (!tw_initialized) //GUI has to be initialized before any displaying occurs (on_create calls post_redisplay by default)
  {
    if (TwInit(TW_OPENGL, NULL))
    {
      tw_initialized = true;
      create_setup_bar();
    }
    else {
      debug_log("E TW init failed: %s\n", TwGetLastError());
    }
  }
# endif

  View::on_create();
}

void ScalarView::on_close()
{
  //GUI cleanup
# ifdef ENABLE_VIEWER_GUI
  if (tw_initialized)
  {
    TwDeleteAllBars();
    TwTerminate();
    tw_initialized = false;
  }
# endif

  //OpenGL clenaup
  if (pointed_node_widget != 0)
  {
    glDeleteLists(pointed_node_widget, 1);
    pointed_node_widget = 0;
  }
  if (selected_node_widget != 0)
  {
    glDeleteLists(selected_node_widget, 1);
    selected_node_widget = 0;
  }
  if (element_id_widget != 0)
  {
    glDeleteLists(element_id_widget, 1);
    element_id_widget = 0;
  }

  //call of parent implementation
  View::on_close();
}

void ScalarView::create_setup_bar()
{
#ifdef ENABLE_VIEWER_GUI
  TwBar* tw_bar = TwNewBar("View setup");

  //contours
  TwAddVarRW(tw_bar, "contours", TW_TYPE_BOOLCPP, &contours, " group=Contour2D label='Show contours'");
  TwAddVarRW(tw_bar, "cont_orig", TW_TYPE_DOUBLE, &cont_orig, " group=Contour2D label='Begin'");
  TwAddVarRW(tw_bar, "cont_step", TW_TYPE_DOUBLE, &cont_step, " group=Contour2D label='Step'");
  TwAddVarRW(tw_bar, "cont_color", TW_TYPE_COLOR3F, &cont_color, " group=Contour2D label='Color'");

  //mesh
  TwAddVarRW(tw_bar, "show_values", TW_TYPE_BOOLCPP, &show_values, " group=Elements2D label='Show Value'");
  TwAddVarRW(tw_bar, "show_edges", TW_TYPE_BOOLCPP, &show_edges, " group=Elements2D label='Show Edges'");
  TwAddVarRW(tw_bar, "show_element_info", TW_TYPE_BOOLCPP, &show_element_info, " group=Elements2D label='Show ID'");
  TwAddVarRW(tw_bar, "edges_color", TW_TYPE_COLOR3F, &edges_color, " group=Elements2D label='Edge color'");

  tw_setup_bar = tw_bar;
#endif
}

void ScalarView::show(MeshFunction* sln, double eps, int item,
                      MeshFunction* xdisp, MeshFunction* ydisp, double dmult)
{
  lin.lock_data();

  double max_abs = range_auto ? -1.0 : std::max(fabs(range_min), fabs(range_max));
  lin.process_solution(sln, item, eps, max_abs, xdisp, ydisp, dmult);
  
  //initialize mesh nodes for displaying and selection
  init_vertex_nodes(sln->get_mesh());

  //initialize element info
  init_element_info(sln->get_mesh());

  if (mode3d) calculate_normals();
  if (range_auto) { range_min = lin.get_min_value();
                    range_max = lin.get_max_value(); }
  glut_init();
  update_layout();

  if (window_id < 0)
  {
    center_mesh(lin.get_vertices(), lin.get_num_vertices());
    center_3d_mesh();
    reset_3d_view();
  }

  create();

  lin.unlock_data();
  verbose("ScalarView::show(): min=%g max=%g", lin.get_min_value(), lin.get_max_value());
  wait_for_draw();
}

bool ScalarView::compare_vertex_nodes_x(const VertexNodeInfo& a, const VertexNodeInfo& b)
{
  return a.x < b.x;
}

void ScalarView::init_vertex_nodes(Mesh* mesh)
{
  //clear all selections
  pointed_vertex_node = NULL;
#ifdef ENABLE_VIEWER_GUI
  if (tw_initialized)
  {
    vector<VertexNodeInfo>::iterator iter = vertex_nodes.begin();
    while (iter != vertex_nodes.end())
    {
      if (iter->tw_bar != NULL)
        TwDeleteBar((TwBar*)iter->tw_bar);
      iter++;
    }
  }
#endif
  vertex_nodes.clear();

  //count a number of active nodes
  int active_nodes = 0;
  int max_node_id = mesh->get_max_node_id();
  for(int i = 0; i < max_node_id; i++)
  {
    Node* node = mesh->get_node(i);
    if (node->type == 0 && node->used) //inspect only vertex nodes
      active_nodes++;
  }

  //allocate
  vertex_nodes.clear();
  vertex_nodes.reserve(active_nodes);

  //copy
  for(int i = 0; i < max_node_id; i++)
  {
    Node* node = mesh->get_node(i);
    if (node->type == 0 && node->used)
      vertex_nodes.push_back(VertexNodeInfo(node->id, (float)node->x, (float)node->y));
  }

  //sort
  std::sort(vertex_nodes.begin(), vertex_nodes.end(), compare_vertex_nodes_x);
}

ScalarView::VertexNodeInfo* ScalarView::find_nearest_node_in_range(float x, float y, float radius)
{
  VertexNodeInfo node_info(-1, x - radius, y); //right side of the widget
  std::vector<VertexNodeInfo>::iterator found_iter = std::lower_bound(vertex_nodes.begin(), vertex_nodes.end(), node_info, compare_vertex_nodes_x);
  std::vector<VertexNodeInfo>::iterator found_nearest = vertex_nodes.end();
  float found_nearest_dist = -1;
  while (found_iter != vertex_nodes.end() && abs(found_iter->x - x) <= radius)
  {
    if (abs(found_iter->y - y) <= radius)
    {
      float dist = min(found_iter->x - x, found_iter->y - y);
      if (found_nearest_dist < 0 || dist < found_nearest_dist)
      {
        found_nearest_dist = dist;
        found_nearest = found_iter;
      }
    }
    found_iter++;
  }

  if (found_nearest != vertex_nodes.end())
    return found_nearest.operator->();
  else
    return NULL;
}

void ScalarView::draw_single_vertex_node(const VertexNodeInfo& node)
{
  //prepare environment
  glPushMatrix();
  glTranslatef(node.x, node.y, 0.0f);
  glScalef(1/(float)scale, 1/(float)scale, 1.0f);

  //prepare number
  void* font = GLUT_BITMAP_HELVETICA_10;
  unsigned char buffer[128];
  sprintf((char*)buffer, "%d", node.id);
  int width_px = glutBitmapLength(font, buffer);
  int height_px = glutBitmapHeight(font);

  //draw target
  if (width_px > (2*node_pixel_radius))
  {
    glPushMatrix();
    float coef = width_px / (2.0f * node_pixel_radius);
    glScalef(coef, coef, 1.0f);
    glCallList(selected_node_widget);
    glPopMatrix();
  }
  else
    glCallList(selected_node_widget);

  //draw number
  glTranslatef(-width_px/2.0f, -height_px/3.0f, 0.0f);
  glColor3f(1.0f, 1.0f, 1.0f);
  glRasterPos2f(0.0f, 0.0f);
  glutBitmapString(font, buffer);

  //clear environment
  glPopMatrix();
}

void ScalarView::draw_vertex_nodes()
{
  //create widgets
  create_nodes_widgets();

  //prepare environment
  glMatrixMode(GL_MODELVIEW);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_LIGHTING);

  //draw selected nodes
  glDisable(GL_BLEND);
  std::vector<VertexNodeInfo>::const_iterator iter = vertex_nodes.begin();
  while (iter != vertex_nodes.end())
  {
    if (iter->selected)
      draw_single_vertex_node(*iter);
    iter++;
  }

  //draw a node under cursor
  if (pointed_vertex_node != NULL)
  {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPushMatrix();
    glTranslatef(pointed_vertex_node->x, pointed_vertex_node->y, 0.0f);
    glScalef(1/(float)scale, 1/(float)scale, 1.0f);
    glCallList(pointed_node_widget);
    glPopMatrix();
    glDisable(GL_BLEND);
  }
}

void ScalarView::create_nodes_widgets()
{
  //pointed node
  if (pointed_node_widget == 0)
  {
    pointed_node_widget  = glGenLists(1);
    glNewList(pointed_node_widget, GL_COMPILE);
    {
      //background
      glColor4f(1.0f, 1.0f, 1.0f, 0.5f);
      glBegin(GL_TRIANGLE_FAN);
      glVertex2f(0.0f, 0.0f);
      float radius = 1.3f * node_pixel_radius;
      for(int i = 0; i < node_widget_vert_cnt; i++)
      {
        float angle = (float)((i / (float)(node_widget_vert_cnt-1)) * (2.0f * M_PI));
        glVertex2f(radius * cos(angle), radius * sin(angle));
      }
      glEnd();

      //foreground
      glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
      glBegin(GL_LINE_STRIP);
      radius = (float)node_pixel_radius;
      for(int i = 0; i < node_widget_vert_cnt; i++)
      {
        float angle = (float)((i / (float)(node_widget_vert_cnt-1)) * (2.0f * M_PI));
        glVertex2f(radius * cos(angle), radius * sin(angle));
      }
      glEnd();
    }
    glEndList();
  }

  //selected mesh node
  if (selected_node_widget == 0)
  {
    selected_node_widget = glGenLists(1);
    glNewList(selected_node_widget, GL_COMPILE);
    {
      //background
      glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
      glBegin(GL_TRIANGLE_FAN);
      glVertex2f(0.0f, 0.0f);
      float radius = 1.3f * node_pixel_radius;
      for(int i = 0; i < node_widget_vert_cnt; i++)
      {
        float angle = (float)((i / (float)(node_widget_vert_cnt-1)) * (2.0f * M_PI));
        glVertex2f(radius * cos(angle), radius * sin(angle));
      }
      glEnd();

      //foreground
      glColor4f(0.4f, 0.4f, 0.2f, 1.0f);
      glBegin(GL_TRIANGLE_FAN);
      glVertex2f(0.0f, 0.0f);
      radius = (float)node_pixel_radius;
      for(int i = 0; i < node_widget_vert_cnt; i++)
      {
        float angle = (float)((i / (float)(node_widget_vert_cnt-1)) * (2.0f * M_PI));
        glVertex2f(radius * cos(angle), radius * sin(angle));
      }
      glEnd();
    }
    glEndList();
  }
}

void ScalarView::init_element_info(Mesh* mesh)
{
  //cleanup
  element_infos.clear();
  
  //check how many active elements are neccessary and prepare space for them
  element_infos.reserve(mesh->get_num_active_elements());

  //build element info
  Element *element = NULL;
  for_all_active_elements(element, mesh)
  {
    double sum_x = 0.0, sum_y = 0.0;
    double max_x, max_y, min_x, min_y;
    max_x = min_x = element->vn[0]->x;
    max_y = min_y = element->vn[0]->y;
    for(int i = 0; i < element->nvert; i++)
    {
      sum_x += element->vn[i]->x;
      sum_y += element->vn[i]->y;

      if (element->vn[i]->x > max_x)
        max_x = element->vn[i]->x;
      if (element->vn[i]->x < min_x)
        min_x = element->vn[i]->x;
      if (element->vn[i]->y > max_y)
        max_y = element->vn[i]->y;
      if (element->vn[i]->y < min_y)
        min_y = element->vn[i]->y;
    }
    element_infos.push_back(ElementInfo(element->id,
      (float)(sum_x / element->nvert), (float)(sum_y / element->nvert),
      (float)(max_x - min_x), (float)(max_y - min_y)));
  }
}

void ScalarView::draw_element_infos_2d()
{
  //create widgets
  create_element_info_widgets();

  //prepare environment
  glMatrixMode(GL_MODELVIEW);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_LIGHTING);
  glDisable(GL_BLEND);

  //draw element IDs
  vector<ElementInfo>::const_iterator iter = element_infos.begin();
  while (iter != element_infos.end())
  {
    //check element dimension in pixels
    float width = (float)(iter->width * scale);
    float height = (float)(iter->height * scale);

    //draw if AABB of element is large enough
    if (width > 3*node_pixel_radius && height > 3*node_pixel_radius)
    {
      //prepare environment
      glPushMatrix();
      glTranslatef(iter->x, iter->y, 0.0f);
      glScalef(1/(float)scale, 1/(float)scale, 1.0f);

      //prepare number
      void* font = GLUT_BITMAP_HELVETICA_10;
      unsigned char buffer[128];
      sprintf((char*)buffer, "%d", iter->id);
      int width_px = glutBitmapLength(font, buffer);
      int height_px = glutBitmapHeight(font);

      //draw background
      glCallList(element_id_widget);

      //draw number
      glTranslatef(-width_px/2.0f, -height_px/3.0f, 0.0f);
      glColor3f(1.0f, 1.0f, 1.0f);
      glRasterPos2f(0.0f, 0.0f);
      glutBitmapString(font, buffer);

      //clear environment
      glPopMatrix();
    }

    //advance to next
    iter++;
  }
}

void ScalarView::create_element_info_widgets()
{
  if (element_id_widget == 0)
  {
    element_id_widget = glGenLists(1);
    glNewList(element_id_widget, GL_COMPILE);
    {
      glBegin(GL_QUADS);

      //background
      float radius = 1.3f * node_pixel_radius;
      glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
      glVertex2f(-radius, radius);
      glVertex2f( radius, radius);
      glVertex2f( radius,-radius);
      glVertex2f(-radius,-radius);

      //foreground
      radius = (float)node_pixel_radius;
      glColor4f(0.2f, 0.2f, 0.4f, 1.0f);
      glVertex2f(-radius, radius);
      glVertex2f( radius, radius);
      glVertex2f( radius,-radius);
      glVertex2f(-radius,-radius);

      glEnd();
    }
    glEndList();
  }
}

void ScalarView::reset_3d_view()
{
  xrot = 40.0;
  yrot = xtrans = ytrans = 0.0;
  ztrans = -3.0;
  yscale = xzscale;
}


void ScalarView::show_contours(double step, double orig)
{
  if (step == 0.0) error("'step' cannot be zero.");
  contours = true;
  cont_orig = orig;
  cont_step = step;
  set_palette_filter(true);
  post_redisplay();
}


static double my_ceil(double x)
{
  double y = ceil(x);
  if (y > x) return y;
  return y + 1.0;
}


void ScalarView::draw_tri_contours(double3* vert, int3* tri)
{
  // sort the vertices by their value, keep track of the permutation sign
  int i, idx[3], perm = 0;
  memcpy(idx, tri, sizeof(idx));
  for (i = 0; i < 2; i++)
  {
    if (vert[idx[0]][2] > vert[idx[1]][2]) { std::swap(idx[0], idx[1]); perm++; }
    if (vert[idx[1]][2] > vert[idx[2]][2]) { std::swap(idx[1], idx[2]); perm++; }
  }
  if (fabs(vert[idx[0]][2] - vert[idx[2]][2]) < 1e-3 * fabs(cont_step)) return;

  // get the first (lowest) contour value
  double val = vert[idx[0]][2];
  val = my_ceil((val - cont_orig) / cont_step) * cont_step + cont_orig;

  int l1 = 0, l2 = 1;
  int r1 = 0, r2 = 2;
  while (val < vert[idx[r2]][2])
  {
    double ld = vert[idx[l2]][2] - vert[idx[l1]][2];
    double rd = vert[idx[r2]][2] - vert[idx[r1]][2];

    // draw a slice of the triangle
    while (val < vert[idx[l2]][2])
    {
      double lt = (val - vert[idx[l1]][2]) / ld;
      double rt = (val - vert[idx[r1]][2]) / rd;
      
      double x1 = (1.0 - lt) * vert[idx[l1]][0] + lt * vert[idx[l2]][0];
      double y1 = (1.0 - lt) * vert[idx[l1]][1] + lt * vert[idx[l2]][1]; 
      double x2 = (1.0 - rt) * vert[idx[r1]][0] + rt * vert[idx[r2]][0];
      double y2 = (1.0 - rt) * vert[idx[r1]][1] + rt * vert[idx[r2]][1]; 
      
      if (perm & 1) { glVertex2d(x1, y1); glVertex2d(x2, y2); }
               else { glVertex2d(x2, y2); glVertex2d(x1, y1); }

      val += cont_step;
    }
    l1 = 1;
    l2 = 2;
  }
}


void ScalarView::on_display()
{
  int i, j;

  lin.lock_data();
  double3* vert = lin.get_vertices();
  int3* tris = lin.get_triangles();
  int3* edges = lin.get_edges();

  // value range
  double min = range_min, max = range_max;
  if (range_auto) { min = lin.get_min_value(); max = lin.get_max_value(); }
  double irange = 1.0 / (max - min);
  // special case: constant solution
  if (fabs(min - max) < 1e-8) { irange = 1.0; min -= 0.5; }

  glPolygonMode(GL_FRONT_AND_BACK, pmode ? GL_LINE : GL_FILL);

  if (!mode3d)
  {
    set_ortho_projection();
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    
    // setup transformation (follows View::transform_x and View::transform_y)
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslated(center_x, center_y, 0.0);
    glScaled(1.0, -1.0, 1.0);
    glTranslated(trans_x, trans_y, 0.0);
    glScaled(scale, scale, 1.0);
    
    // draw all triangles
    if (show_values)
    {
      if (gl_pallete_tex_id == 0)
        debug_log("E palette texture ID is zero, palette not set\n");
      glEnable(GL_TEXTURE_1D);
      glBindTexture(GL_TEXTURE_1D, gl_pallete_tex_id);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
      glColor3f(0, 1, 0);
      glBegin(GL_TRIANGLES);
      for (i = 0; i < lin.get_num_triangles(); i++)
      {
        if (finite(vert[tris[i][0]][2]) && finite(vert[tris[i][1]][2]) && finite(vert[tris[i][2]][2]))
        {
          glTexCoord2d((vert[tris[i][0]][2] - min) * irange * tex_scale + tex_shift, 0.0);
          glVertex2d(vert[tris[i][0]][0], vert[tris[i][0]][1]);
          
          glTexCoord2d((vert[tris[i][1]][2] - min) * irange * tex_scale + tex_shift, 0.0);
          glVertex2d(vert[tris[i][1]][0], vert[tris[i][1]][1]);
          
          glTexCoord2d((vert[tris[i][2]][2] - min) * irange * tex_scale + tex_shift, 0.0);
          glVertex2d(vert[tris[i][2]][0], vert[tris[i][2]][1]);
        }
      }
      glEnd();
    }
    
    // draw contours
    glDisable(GL_TEXTURE_1D);
    if (contours)
    {
      glColor3fv(cont_color);
      glBegin(GL_LINES);
      for (i = 0; i < lin.get_num_triangles(); i++)
      {
        if (finite(vert[tris[i][0]][2]) && finite(vert[tris[i][1]][2]) && finite(vert[tris[i][2]][2]))
        {
          draw_tri_contours(vert, &tris[i]);
        }
      }
      glEnd();
    }

    // draw all edges
    if (show_edges)
    {
      glColor3fv(edges_color);
      glBegin(GL_LINES);
      for (i = 0; i < lin.get_num_edges(); i++) // todo: we could draw only left-right, top-bottom ones
      {
        if (lines || edges[i][2] != 0)
        {
          glVertex2d(vert[edges[i][0]][0], vert[edges[i][0]][1]);
          glVertex2d(vert[edges[i][1]][0], vert[edges[i][1]][1]);
        }
      }
      glEnd();
    }

    //draw element IDS
    if (show_element_info)
      draw_element_infos_2d();

    //draw nodes and node info
    if (show_edges)
      draw_vertex_nodes();

    //cleanup
    glPopMatrix();
  }
  else
  {
    set_3d_projection(50, 0.05, 10.0);
    glClear(GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    glLoadIdentity();
    glTranslated(xtrans, ytrans, ztrans);
    double xr = xrot;
    glRotated(xr, 1, 0, 0);
    glRotated(yrot, 0, 1, 0);

    // draw the surface
    glEnable(GL_LIGHTING);
    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D, gl_pallete_tex_id);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);
    glBegin(GL_TRIANGLES);
    for (i = 0; i < lin.get_num_triangles(); i++)
    {
      for (j = 0; j < 3; j++)
      {
        glNormal3d(normals[tris[i][j]][0], normals[tris[i][j]][2], -normals[tris[i][j]][1]);
        glTexCoord2d((vert[tris[i][j]][2] - min) * irange * tex_scale + tex_shift, 0.0);
        glVertex3d((vert[tris[i][j]][0] - xctr) * xzscale,
                   (vert[tris[i][j]][2] - yctr) * yscale,
                  -(vert[tris[i][j]][1] - zctr) * xzscale);
      }
    }
    glEnd();
    glDisable(GL_POLYGON_OFFSET_FILL);

    // draw edges
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_1D);
    if (lines)
    {
      glColor3f(0.5, 0.5, 0.5);
      //glColor3f(0.0, 0.0, 0.0);
      glBegin(GL_LINES);
      for (i = 0; i < lin.get_num_edges(); i++)
      {
        glVertex3d((vert[edges[i][0]][0] - xctr) * xzscale,
                   (vert[edges[i][0]][2] - yctr) * yscale,
                  -(vert[edges[i][0]][1] - zctr) * xzscale);
        glVertex3d((vert[edges[i][1]][0] - xctr) * xzscale,
                   (vert[edges[i][1]][2] - yctr) * yscale,
                  -(vert[edges[i][1]][1] - zctr) * xzscale);
      }
      glEnd();
    }

    // draw boundary edges
    glColor3f(0, 0, 0);
    glBegin(GL_LINES);
    for (i = 0; i < lin.get_num_edges(); i++)
    {
      if (edges[i][2])
      {
        glVertex3d((vert[edges[i][0]][0] - xctr) * xzscale, -yctr * yscale,
                  -(vert[edges[i][0]][1] - zctr) * xzscale);
        glVertex3d((vert[edges[i][1]][0] - xctr) * xzscale, -yctr * yscale,
                  -(vert[edges[i][1]][1] - zctr) * xzscale);
      }
    }
    glEnd();
  }

  //draw TW
#ifdef ENABLE_VIEWER_GUI
  if (tw_initialized)
    TwDraw();
#endif
  
  lin.unlock_data();
}

static inline void normalize(double& x, double& y, double& z)
{
  double l = 1.0 / sqrt(sqr(x) + sqr(y) + sqr(z));
  x *= l; y *= l; z *= l;
}


void ScalarView::calculate_normals()
{
  lin.lock_data();
  int nv = lin.get_num_vertices();
  int nt = lin.get_num_triangles();
  double3* vert = lin.get_vertices();
  int3* tris = lin.get_triangles();

  if (normals != NULL) delete [] normals;
  normals = new double3[nv];
  memset(normals, 0, nv * sizeof(double3));

  for (int i = 0; i < nt; i++)
  {
    double ax = (vert[tris[i][1]][0] - vert[tris[i][0]][0]) * xzscale;
    double ay = (vert[tris[i][1]][1] - vert[tris[i][0]][1]) * xzscale;
    double az = (vert[tris[i][1]][2] - vert[tris[i][0]][2]) * yscale;

    double bx = (vert[tris[i][2]][0] - vert[tris[i][0]][0]) * xzscale;
    double by = (vert[tris[i][2]][1] - vert[tris[i][0]][1]) * xzscale;
    double bz = (vert[tris[i][2]][2] - vert[tris[i][0]][2]) * yscale;

    double nx = ay * bz - az * by;
    double ny = az * bx - ax * bz;
    double nz = ax * by - ay * bx;
    normalize(nx, ny, nz);

    for (int j = 0; j < 3; j++)
    {
      normals[tris[i][j]][0] += nx;
      normals[tris[i][j]][1] += ny;
      normals[tris[i][j]][2] += nz;
    }
  }

  for (int i = 0; i < nv; i++)
    normalize(normals[i][0], normals[i][1], normals[i][2]);

  lin.unlock_data();
}


void ScalarView::center_3d_mesh()
{
  lin.lock_data();
  double3* vert = lin.get_vertices();

  double xmin = 1e100, xmax = -1e100;
  double zmin = 1e100, zmax = -1e100;
  double ysum = 0.0;
  for (int i = 0; i < lin.get_num_vertices(); i++)
  {
    if (vert[i][0] < xmin) xmin = vert[i][0];
    if (vert[i][0] > xmax) xmax = vert[i][0];
    if (vert[i][1] < zmin) zmin = vert[i][1];
    if (vert[i][1] > zmax) zmax = vert[i][1];
    ysum += vert[i][2];
  }

  xzscale = 1.8 / std::max(xmax - xmin, zmax - zmin);
  xctr = (xmax + xmin) / 2.0;
  zctr = (zmax + zmin) / 2.0;
  yscale  = xzscale;
  yctr = ysum / lin.get_num_vertices();

  lin.unlock_data();
}


void ScalarView::init_lighting()
{
  float light_specular[] = {  1.0f, 1.0f, 1.0f, 1.0f };
  float light_ambient[]  = {  0.3f, 0.3f, 0.3f, 1.0f };
  float light_diffuse[]  = {  1.0f, 1.0f, 1.0f, 1.0f };
  float light_position[] = { -0.5f, 0.5f, 1.0f, 0.0f };
  
  glEnable(GL_LIGHT0);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);

  float material_ambient[]  = { 0.2f, 0.2f, 0.2f, 1.0f };
  float material_diffuse[]  = { 0.8f, 0.8f, 0.8f, 1.0f };
  float material_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };

  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_diffuse);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128);
  glDisable(GL_COLOR_MATERIAL);

  glShadeModel(GL_SMOOTH);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
#if defined(GL_LIGHT_MODEL_COLOR_CONTROL) && defined(GL_SEPARATE_SPECULAR_COLOR)
  glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);
#endif
}


void ScalarView::on_key_down(unsigned char key, int x, int y)
{
  VIEWER_GUI_CALLBACK(TwEventKeyboardGLUT(key, x, y))
  {
    switch (key)
    {
      case 'm':
        lines = !lines;
        post_redisplay();
        break;

      case 'l':
        pmode = !pmode;
        post_redisplay();
        break;

      case 'c':
      {
        if (!mode3d) {
          lin.lock_data();
          center_mesh(lin.get_vertices(), lin.get_num_vertices());
          lin.unlock_data();
        }
        else
        {
          center_3d_mesh();
          reset_3d_view();
        }
        post_redisplay();
        //reset_zoom();
        break;
      }

      case 'f':
        set_palette_filter(pal_filter != GL_LINEAR);
        break;
      
      case 'k':
        contours = !contours;
        post_redisplay();
        break;

      case 'i':
        show_element_info = !show_element_info;
        if (!mode3d)
          post_redisplay();
        break;
      
      case '3':
      {
        mode3d = !mode3d;
        dragging = scaling = false;
        if (mode3d)
        {
          init_lighting();
          if (normals == NULL) calculate_normals();
        }
        post_redisplay();
        break;
      }
      
      case '*':
      case '/':
      {
        if (mode3d)
        {
          if (key == '*') yscale *= 1.1; else yscale /= 1.1;
          calculate_normals();
        }
        else if (contours)
        {
          if (key == '*') cont_step /= 1.1; else cont_step *= 1.1;
        }
        post_redisplay();
        break;
      }

      default:
        View::on_key_down(key, x, y);
        break;
    }
  }
}

void ScalarView::on_special_key(int key, int x, int y)
{
  VIEWER_GUI_CALLBACK(TwEventSpecialGLUT(key, x, y))
  {
    View::on_special_key(key, x, y);
  }
}


void ScalarView::on_mouse_move(int x, int y)
{
  VIEWER_GUI_CALLBACK(TwEventMouseMotionGLUT(x, y))
  {
    if (dragging)
    {
      yrot += 0.4 * (x - mouse_x);
      xrot += 0.4 * (y - mouse_y);

      if (xrot < -90) xrot = -90;
      else if (xrot > 90) xrot = 90;
    }
    else if (scaling)
    {
      if (dragging)
      {
        yrot += 0.4 * (x - mouse_x);
        xrot += 0.4 * (y - mouse_y);
        
        if (xrot < -90) xrot = -90;
        else if (xrot > 90) xrot = 90;
      }
      else if (scaling)
      {
        ztrans += 0.01 * (mouse_y - y);
        if (ztrans > -0.25) ztrans = -0.25;
        else if (ztrans < -7) ztrans = -7;
      }
      else
      {
        xtrans += 0.002 * (x - mouse_x);
        ytrans += 0.002 * (mouse_y - y);
      }
        
      post_redisplay();
      mouse_x = x;
      mouse_y = y;
      return;
    }
    else if (!mode3d && show_edges && !dragging && !scaling && !panning)
    {
      VertexNodeInfo* found_node = find_nearest_node_in_range((float)untransform_x(x), (float)untransform_y(y), (float)(node_pixel_radius / scale));
      if (found_node != pointed_vertex_node)
        pointed_vertex_node = found_node;
    }

    post_redisplay();
  }
}

void ScalarView::on_left_mouse_down(int x, int y)
{
  VIEWER_GUI_CALLBACK(TwEventMouseButtonGLUT(GLUT_LEFT_BUTTON, GLUT_DOWN, x, y))
  {
    View::on_left_mouse_down(x, y);
  }
}

void ScalarView::on_left_mouse_up(int x, int y)
{
  VIEWER_GUI_CALLBACK(TwEventMouseButtonGLUT(GLUT_LEFT_BUTTON, GLUT_UP, x, y))
  {
    View::on_left_mouse_up(x, y);
  }
}

void ScalarView::on_right_mouse_down(int x, int y)
{
  VIEWER_GUI_CALLBACK(TwEventMouseButtonGLUT(GLUT_RIGHT_BUTTON, GLUT_DOWN, x, y))
  {
    //handle node selection
    if (pointed_vertex_node != NULL)
    {
      if (pointed_vertex_node->selected) //deselct node
      {
# ifdef ENABLE_VIEWER_GUI
        if (pointed_vertex_node->tw_bar != NULL)
        {
          TwDeleteBar((TwBar*)pointed_vertex_node->tw_bar);
          pointed_vertex_node->tw_bar = NULL;
        }
# endif
        pointed_vertex_node->selected = false;
      }
      else { //select node
# ifdef ENABLE_VIEWER_GUI
        char buffer[1024];
        sprintf(buffer, "Vertex node: %d", pointed_vertex_node->id);
        TwBar* bar = TwNewBar(buffer);
        TwAddVarRO(bar, "X", TW_TYPE_FLOAT, &pointed_vertex_node->x, NULL);
        TwAddVarRO(bar, "Y", TW_TYPE_FLOAT, &pointed_vertex_node->y, NULL);
        int intTmp = 0; TwSetParam(bar, NULL, "iconifiable", TW_PARAM_INT32, 1, &intTmp);
        pointed_vertex_node->tw_bar = bar;
# endif
        pointed_vertex_node->selected = true;
      }
      post_redisplay();
    }
    else { //avoid mixing of functionality
      View::on_right_mouse_down(x, y);
    }
  }
}

void ScalarView::on_right_mouse_up(int x, int y)
{
  VIEWER_GUI_CALLBACK(TwEventMouseButtonGLUT(GLUT_RIGHT_BUTTON, GLUT_UP, x, y))
  {
    View::on_right_mouse_up(x, y);
  }
}

void ScalarView::on_middle_mouse_down(int x, int y)
{
  VIEWER_GUI_CALLBACK(TwEventMouseButtonGLUT(GLUT_MIDDLE_BUTTON, GLUT_DOWN, x, y))
  {
    if (!mode3d) return;
    dragging = scaling = false;
    panning = true;
  }
}

void ScalarView::on_middle_mouse_up(int x, int y)
{
  VIEWER_GUI_CALLBACK(TwEventMouseButtonGLUT(GLUT_MIDDLE_BUTTON, GLUT_UP, x, y))
  {
    panning = false;
  }
}

void ScalarView::on_reshape(int width, int height)
{
  VIEWER_GUI(TwWindowSize(width, height));

  View::on_reshape(width, height);
}


//// load & save ///////////////////////////////////////////////////////////////////////////////////

void ScalarView::load_data(const char* filename)
{
  lin.lock_data();
  lin.load_data(filename);
  if (mode3d) calculate_normals();
  if (range_auto) { range_min = lin.get_min_value();
                    range_max = lin.get_max_value(); }
  glut_init();
  update_layout();

  if (window_id < 0)
  {
    center_mesh(lin.get_vertices(), lin.get_num_vertices());
    center_3d_mesh();
    reset_3d_view();
  }

  create();
  lin.unlock_data();
  wait_for_draw();
}


void ScalarView::save_data(const char* filename)
{
  if (lin.get_num_triangles() <= 0)
    error("No data to save.");
  lin.save_data(filename);
}


void ScalarView::save_numbered(const char* format, int number)
{
  char buffer[1000];
  sprintf(buffer, format, number);
  save_data(buffer);
}


const char* ScalarView::get_help_text() const
{
  return
  "ScalarView\n\n"
  "Controls:\n"
  "  Left mouse - pan\n"
  "  Right mouse - zoom\n"
  "  3 - toggle 3D mode\n"
  "  C - center image\n"
  "  F - toggle smooth palette\n"
  "  H - render high-quality frame\n"
  "  K - toggle contours\n"
  "  M - toggle mesh\n"
  "  I - toggle element IDs (2d only)\n"
  "  P - cycle palettes\n"
  "  S - save screenshot\n"
  "  * - increase contour density\n"
  "  / - decrease contour density\n"
  "  F1 - this help\n"
  "  Esc, Q - quit\n\n"
  "3D mode:\n"
  "  Left mouse - rotate\n"
  "  Right mouse - zoom\n"
  "  Middle mouse - pan\n"
  "  * - increase Z scale\n"
  "  / - decrease Z scale";
}


#endif // NOGLUT
