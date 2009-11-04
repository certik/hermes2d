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

#ifndef __HERMES2D_ORDER_VIEW_H
#define __HERMES2D_ORDER_VIEW_H

#include "view.h"

// you can define NOGLUT to turn off all OpenGL stuff in Hermes2D
#ifndef NOGLUT

/// \brief Displays the polynomial degrees of elements.
///
/// OrderView is a tool for displaying the polynomial degrees of the elements in a space.
///
class PUBLIC_API OrderView : public View
{
public:

  OrderView(const char* title = "OrderView", DEFAULT_WINDOW_POS);

  void show(Space* space);

  void load_data(const char* filename);
  void save_data(const char* filename);
  void save_numbered(const char* format, int number);

protected:

  Orderizer ord;
  bool b_orders;

  int num_boxes, order_min;
  const char* box_names[11];
  char text_buffer[500];
  float order_colors[11][3];

  void init_order_palette();

  virtual void on_display();
  virtual void on_key_down(unsigned char key, int x, int y);
  virtual void scale_dispatch();
  virtual int measure_scale_labels();
  virtual const char* get_help_text() const;

};

#else // NOGLUT

class PUBLIC_API OrderView : public View
{
public:
  OrderView(const char* title = "OrderView", DEFAULT_WINDOW_POS) {}
  void show(Space* space)
     { info("OrderView: Hermes2D compiled without OpenGL support, skipping visualization."); }
  void load_data(const char* filename) {}
  void save_data(const char* filename) {}
  void save_numbered(const char* format, int number) {}
};

#endif // NOGLUT

#endif
