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

// $Id: view.h 1086 2008-10-21 09:05:44Z jakub $

#ifndef __HERMES2D_VIEW_H
#define __HERMES2D_VIEW_H

#include "../common.h"
#include "../linear.h"

// default window size and position
#define DEFAULT_WINDOW_POS  int x = -1, int y = -1, int width = 1000, int height = 800

// you can define NOGLUT to turn off all OpenGL stuff in Hermes2D
#ifndef NOGLUT

/// View palette type.
enum ViewPaletteType {
  H2DV_PT_DEFAULT = -1, ///< Default palette. Depends on viewer.
  H2DV_PT_HUESCALE = 0, ///< A palette based on hue scale.
  H2DV_PT_GRAYSCALE = 1, ///< Greyscale.
  H2DV_PT_INVGRAYSCALE = 2, ///< Inverted grayscale.
  H2DV_PT_MAX_ID = 3 ///< Maximum ID of view palette type.
};

/// \brief Represents a simple visualization window.
///
/// View is a base class providing a simple OpenGL visualization window.
/// Its task is to define basic functionality, such as the ability of the
/// window to be responsive even when the main program thread is busy
/// with calculations (ie., the windows are run in a background thread),
/// to provide zooming and panning capabilities for use by the descendant
/// classes, etc.
///
class PUBLIC_API View
{
public:

  View(const char* title, int x, int y, int width, int height);
  virtual ~View();

  int  create();
  void close();
  void refresh(); ///< Refreshes views

  /// Changes the window name (in its title-bar) to 'title'.
  void set_title(const char* title);

  void set_min_max_range(double min, double max);
  void auto_min_max_range();
  void get_min_max_range(double& min, double& max);

  void show_scale(bool show = true);
  void set_scale_position(int horz, int vert);
  void set_scale_size(int width, int height, int numticks);
  void set_scale_format(const char* fmt);
  void fix_scale_width(int width = 80);

  /// Saves the current content of the window to a .BMP file.
  /// If 'high_quality' is true, an anti-aliased frame is rendered and saved.
  void save_screenshot(const char* bmpname, bool high_quality = false);
  /// Like save_screenshot(), but forms the file name in printf-style using the 'number'
  /// parameter, e.g., format="screen%03d.bmp" and number=5 gives the file name "screen005.bmp".
  void save_numbered_screenshot(const char* format, int number, bool high_quality = false);

  void set_palette(ViewPaletteType type);
  void set_num_palette_steps(int num);
  void set_palette_filter(bool linear);

  void wait_for_keypress();
  void wait_for_keypress(const char* text);
  void wait_for_close();
  void wait_for_draw();

  static void wait(const char* text = NULL); ///< Closes all views at once.

protected: //FPS measurement
#define FPS_FRAME_SIZE 5
  double rendering_frames[FPS_FRAME_SIZE]; ///< time spend in rendering of frames [in ms]
  int rendering_frames_top; ///< the new location of the next FPS
  void draw_fps(); ///< draws current FPS
  static double get_tick_count(); ///< returns a current time [in ms]

protected:

  virtual void clear_background();

  virtual void on_create(int output_id);
  virtual void on_display() = 0;
  virtual void on_reshape(int width, int height);
  virtual void on_mouse_move(int x, int y);
  virtual void on_left_mouse_down(int x, int y);
  virtual void on_left_mouse_up(int x, int y);
  virtual void on_left_mouse_double_click(int x, int y) {}
  virtual void on_right_mouse_down(int x, int y);
  virtual void on_right_mouse_up(int x, int y);
  virtual void on_right_mouse_double_click(int x, int y) {}
  virtual void on_middle_mouse_down(int x, int y) {}
  virtual void on_middle_mouse_up(int x, int y) {}
  virtual void on_middle_mouse_double_click(int x, int y) {}
  virtual void on_key_down(unsigned char key, int x, int y);
  virtual void on_special_key(int key, int x, int y);
  virtual void on_entry(int state) {}
  virtual void on_close();


  template<class TYPE> void center_mesh(TYPE* vertices, int nvert);
  virtual void get_palette_color(double x, float* gl_color); ///< Fills gl_color with palette color. Assumes that gl_color points to a vector of three components (RGB).

protected:

  std::string title;
  int output_id;
  int output_x, output_y, output_width, output_height;
  float jitter_x, jitter_y;
  bool hq_frame, frame_ready;

  double scale, log_scale, trans_x, trans_y;
  double center_x, center_y;
  int margin, lspace, rspace;
  int mouse_x, mouse_y;
  int scx, scy;
  double objx, objy;
  bool dragging, scaling;

  ViewPaletteType pal_type;
  int pal_steps, pal_filter;
  double tex_scale, tex_shift;
  bool range_auto;
  double range_min, range_max;

  bool b_scale, b_help;
  bool scale_focused, scale_dragging;
  int pos_horz, pos_vert;
  int scale_x, scale_y;
  int scale_width, scale_height, labels_width;
  int scale_numticks, scale_box_height, scale_box_skip;
  char scale_fmt[20];
  int scale_fixed_width;

  bool want_screenshot;
  static int screenshot_no;
  std::string screenshot_filename;

  void update_scale();
  void update_log_scale();

protected: //OpenGL
  unsigned int gl_pallete_tex_id;

protected: //internal functions
  double transform_x(double x) { return (x * scale + trans_x) + center_x; }
  double transform_y(double y) { return center_y - (y * scale + trans_y); }
  double untransform_x(double x) { return (x - center_x - trans_x) / scale; }
  double untransform_y(double y) { return (center_y - y - trans_y) / scale; }

  void pre_display();
  void display_antialiased();

  void set_ortho_projection(bool no_jitter = false);
  void set_3d_projection(int fov, double znear, double zfar);

  void draw_text(double x, double y, const char* text, int align = -1);
  int  get_text_width(const char* text);

  char *get_screenshot_file_name();
  void save_screenshot_internal(const char* filename);

  virtual void scale_dispatch();
  virtual int measure_scale_labels();
  void draw_continuous_scale(char* title, bool righttext);
  void draw_discrete_scale(int numboxes, const char* boxnames[], const float boxcolors[][3]);

  void create_gl_palette(); ///< Creates pallete texture in OpenGL. Assumes that view_sync is locked.
  void update_tex_adjust();
  void update_layout();

  void draw_help();
  virtual const char* get_help_text() const { return ""; }

  friend void on_display_stub(void);
  friend void on_reshape_stub(int, int);
  friend void on_mouse_move_stub(int, int);
  friend void on_mouse_click_stub(int, int, int, int);
  friend void on_key_down_stub(unsigned char, int, int);
  friend void on_special_key_stub(int, int, int);
  friend void on_entry_stub(int);
  friend void on_idle_stub();
  friend void on_close_stub();
  friend int add_view_in_thread(void*);
  friend int remove_view_in_thread(void*);
};

// this has to be here, unfortunatelly (because the author is not able to use pointers)
template<class TYPE>
void View::center_mesh(TYPE* vertices, int nvert)
{
  if (nvert <= 0) return;

  // get mesh bounding box
  double xmin = 1e10, xmax = -1e10;
  double ymin = 1e10, ymax = -1e10;
  for (int i = 0; i < nvert; i++)
  {
    if (vertices[i][0] < xmin) xmin = vertices[i][0];
    if (vertices[i][0] > xmax) xmax = vertices[i][0];
    if (vertices[i][1] < ymin) ymin = vertices[i][1];
    if (vertices[i][1] > ymax) ymax = vertices[i][1];
  }
  double mesh_width  = xmax - xmin;
  double mesh_height = ymax - ymin;

  double usable_width = output_width - 2*margin - lspace - rspace;
  double usable_height = output_height - 2*margin;

  // align in the proper direction
  if (usable_width / usable_height < mesh_width / mesh_height)
    scale = usable_width / mesh_width;
  else
    scale = usable_height / mesh_height;

  // center
  trans_x = -scale * (xmin + xmax) / 2;
  trans_y = -scale * (ymin + ymax) / 2;

  update_log_scale();
}


#else // NOGLUT

class PUBLIC_API View
{
public:
  View() {}
  View(const char* title, int x, int y, int width, int height) {}
  ~View() {}
  int  create() { return 0; }
  void close() {}
  void set_title(const char* title) {}
  void set_min_max_range(double min, double max) {}
  void auto_min_max_range() {}
  void get_min_max_range(double& min, double& max) {}
  void show_scale(bool show = true) {}
  void set_scale_position(int horz, int vert) {}
  void set_scale_size(int width, int height, int numticks) {}
  void set_scale_format(const char* fmt) {}
  void fix_scale_width(int width = 80) {}
  void save_screenshot(const char* bmpname, bool high_quality = false) {}
  void save_numbered_screenshot(const char* format, int number, bool high_quality = false) {}
  void set_palette(int type) {}
  void set_num_palette_steps(int num) {}
  void set_palette_filter(bool linear) {}
  void wait_for_keypress() {}
  void wait_for_close() {}
  void wait_for_draw() {}
  static void wait(const char* text = NULL) {}
};

#endif // NOGLUT

#endif
