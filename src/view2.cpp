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

// $Id: view2.cpp 1086 2008-10-21 09:05:44Z jakub $

#ifndef NOGLUT

#include <GL/freeglut.h>
#include "common.h"
#include "view.h"


//// ScalarView ////////////////////////////////////////////////////////////////////////////////////

ScalarView::ScalarView(const char* title, int x, int y, int width, int height)
           : View(title, x, y, width, height)
{
  lines = true;
  pmode = mode3d = false;
  normals = NULL;
  panning = false;
  contours = false;
  cont_orig = 0.0;
  cont_step = 1.0;
}

ScalarView::~ScalarView()
{
  if (normals != NULL) delete [] normals;
}


void ScalarView::show(MeshFunction* sln, double eps, int item,
                      MeshFunction* xdisp, MeshFunction* ydisp, double dmult)
{
  lin.lock_data();

  double max_abs = range_auto ? -1.0 : std::max(fabs(range_min), fabs(range_max));
  lin.process_solution(sln, item, eps, max_abs, xdisp, ydisp, dmult);

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


void ScalarView::draw_tri_contours(double3* vert, double2* tvert, int3* tri)
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

      double x1 = (1.0 - lt) * tvert[idx[l1]][0] + lt * tvert[idx[l2]][0];
      double y1 = (1.0 - lt) * tvert[idx[l1]][1] + lt * tvert[idx[l2]][1];
      double x2 = (1.0 - rt) * tvert[idx[r1]][0] + rt * tvert[idx[r2]][0];
      double y2 = (1.0 - rt) * tvert[idx[r1]][1] + rt * tvert[idx[r2]][1];

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

    // transform all vertices
    int nv = lin.get_num_vertices();
    double2* tvert = new double2[nv];
    for (i = 0; i < nv; i++)
    {
      tvert[i][0] = transform_x(vert[i][0]);
      tvert[i][1] = transform_y(vert[i][1]);
    }

    // draw all triangles
    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D, 1);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    glBegin(GL_TRIANGLES);
    for (i = 0; i < lin.get_num_triangles(); i++)
    {
      if (finite(vert[tris[i][0]][2]) && finite(vert[tris[i][1]][2]) && finite(vert[tris[i][2]][2]))
      {
        glTexCoord2d((vert[tris[i][0]][2] - min) * irange * tex_scale + tex_shift, 0.0);
        glVertex2d(tvert[tris[i][0]][0], tvert[tris[i][0]][1]);

        glTexCoord2d((vert[tris[i][1]][2] - min) * irange * tex_scale + tex_shift, 0.0);
        glVertex2d(tvert[tris[i][1]][0], tvert[tris[i][1]][1]);

        glTexCoord2d((vert[tris[i][2]][2] - min) * irange * tex_scale + tex_shift, 0.0);
        glVertex2d(tvert[tris[i][2]][0], tvert[tris[i][2]][1]);
      }
    }
    glEnd();

    // draw contours
    glDisable(GL_TEXTURE_1D);
    if (contours)
    {
      glColor3f(0.0, 0.0, 0.0);
      glBegin(GL_LINES);
      for (i = 0; i < lin.get_num_triangles(); i++)
      {
        if (finite(vert[tris[i][0]][2]) && finite(vert[tris[i][1]][2]) && finite(vert[tris[i][2]][2]))
        {
          draw_tri_contours(vert, tvert, &tris[i]);
        }
      }
      glEnd();
    }

    // draw all edges
    glColor3f(0.5, 0.5, 0.5);
    glBegin(GL_LINES);
    for (i = 0; i < lin.get_num_edges(); i++) // todo: we could draw only left-right, top-bottom ones
    {
      if (lines || edges[i][2] != 0)
      {
        glVertex2d(tvert[edges[i][0]][0], tvert[edges[i][0]][1]);
        glVertex2d(tvert[edges[i][1]][0], tvert[edges[i][1]][1]);
      }
    }
    glEnd();

    delete [] tvert;
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
    glBindTexture(GL_TEXTURE_1D, 1);
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

  lin.unlock_data();
}


void ScalarView::reset_3d_view()
{
  xrot = 40.0;
  yrot = xtrans = ytrans = 0.0;
  ztrans = -3.0;
  yscale = xzscale;
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
  float light_specular[] = {  1.0, 1.0, 1.0, 1.0 };
  float light_ambient[]  = {  0.3, 0.3, 0.3, 1.0 };
  float light_diffuse[]  = {  1.0, 1.0, 1.0, 1.0 };
  float light_position[] = { -0.5, 0.5, 1.0, 0.0 };

  glEnable(GL_LIGHT0);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);

  float material_ambient[]  = { 0.2, 0.2, 0.2, 1.0 };
  float material_diffuse[]  = { 0.8, 0.8, 0.8, 1.0 };
  float material_specular[] = { 1.0, 1.0, 1.0, 1.0 };

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


void ScalarView::on_mouse_move(int x, int y)
{
  if (mode3d && (dragging || scaling || panning))
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
  View::on_mouse_move(x, y);
}


void ScalarView::on_middle_mouse_down(int x, int y)
{
  if (!mode3d) return;
  dragging = scaling = false;
  panning = true;
}

void ScalarView::on_middle_mouse_up(int x, int y)
{
  panning = false;
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
