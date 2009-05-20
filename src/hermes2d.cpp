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

// $Id: hermes2d.cpp 1086 2008-10-21 09:05:44Z jakub $

#include "common.h"
#include "hermes2d.h"


void hermes2d_initialize(int* argc, char* argv[])
{
  warn("this function is deprecated.");
}


void hermes2d_finalize(bool force_quit)
{
  warn("this function is deprecated.");
  #ifndef NOGLUT
      finish_glut_main_loop(force_quit);
  #endif
  //free_ortho_base();
  // TODO: free matrices in curved.cpp
  // TODO: free bases in H1OrthoHp, ...
}
