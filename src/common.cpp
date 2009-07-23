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

// $Id: common.cpp 1086 2008-10-21 09:05:44Z jakub $

#include "common.h"


bool verbose_mode = true;
bool info_mode = true;
bool warn_integration = true;


void __error_fn(const char* fname, const char* msg, ...)
{
  char text[500];
  sprintf(text, "ERROR: %s: %s\n", fname, msg);
  va_list arglist;
  va_start(arglist, msg);
  vfprintf(stderr, text, arglist);
  va_end(arglist);
  exit(1);
}


void __warn_fn(const char* fname, const char* msg, ...)
{
  char text[500];
  sprintf(text, "WARNING: %s: %s\n", fname, msg);
  va_list arglist;
  va_start(arglist, msg);
  vfprintf(stderr, text, arglist);
  va_end(arglist);
}


void __info_fn(const char* msg, ...)
{
  if (!info_mode) return;
  char text[500];
  sprintf(text, "%s\n", msg);
  va_list arglist;
  va_start(arglist, msg);
  vprintf(text, arglist);
  va_end(arglist);
}


void __verbose_fn(const char* msg, ...)
{
  if (!verbose_mode) return;
  char text[500];
  sprintf(text, "%s\n", msg);
  va_list arglist;
  va_start(arglist, msg);
  vprintf(text, arglist);
  va_end(arglist);
}


void __hermes2d_fwrite(const void* ptr, size_t size, size_t nitems, FILE* stream, const char* fname)
{
  if (fwrite(ptr, size, nitems, stream) != nitems || ferror(stream))
    __error_fn(fname, "Error writing to file: %s", strerror(ferror(stream)));
}

void __hermes2d_fread(void* ptr, size_t size, size_t nitems, FILE* stream, const char* fname)
{
  size_t ret = fread(ptr, size, nitems, stream);
  if (ret < nitems)
    __error_fn(fname, "Premature end of file.");
  else if (ferror(stream))
    __error_fn(fname, "Error reading file: %s", strerror(ferror(stream)));
}


//// timing stuff //////////////////////////////////////////////////////////////////////////////////

/*#ifndef WIN32
  #include <sys/time.h>
#else*/
  #include <time.h>
//#endif


static const int max_stack = 50;
static clock_t tick_stack[max_stack];
static int ts_top = 0;


void begin_time() // TODO: make this return wall time on both Linux and Win32
{
  if (ts_top < max_stack) tick_stack[ts_top] = clock();
  else warn("Timing stack overflow.");
  ts_top++;
}


double cur_time()
{
  if (!ts_top) { warn("Called without begin_time()."); return -1.0; }
  if (ts_top >= max_stack) return -1.0;
  return (double) (clock() - tick_stack[ts_top-1]) / CLOCKS_PER_SEC;
}


double end_time()
{
  if (!ts_top) { warn("Called without begin_time()."); return -1.0; }
  ts_top--;
  if (ts_top >= max_stack) return -1.0;
  return (double) (clock() - tick_stack[ts_top]) / CLOCKS_PER_SEC;
}

void throw_exception(char *text)
{
    throw std::runtime_error(text);
}
