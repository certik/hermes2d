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

#ifndef __HERMES2D_VIEWS_OUTPUT_H
#define __HERMES2D_VIEWS_OUTPUT_H

// you can define NOGLUT to turn off all OpenGL stuff in Hermes2D
#ifndef NOGLUT

#include <vector>
class View;

// exported variables
extern pthread_t thread;
extern int num_windows; ///< A number of GLUT windows
extern std::vector<View*> wnd_instance; ///< Instances of viewers. ID is taken from GLUT system.

// exported methods 
extern void glut_init();
extern void finish_glut_main_loop(bool force); //deprecated: don't use
extern int cross_thread_call(int (*function)(void*), void* param); 

#endif // NOGLUT

#endif
