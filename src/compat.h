/*
   This header file contains platform compatibility layer.

   It is included from common.h, so it is automatically included in all hermes
   sources. The implementation of the functions in this file is in the
   src/compat directory.
*/

#ifndef __HERMES2D_FMEMOPEN_H
#define __HERMES2D_FMEMOPEN_H

#include "config.h"

#ifndef HAVE_FMEMOPEN

/* The GNU libc fmemopen
   This function is available by default in GNU libc and on other platforms
   (like OSX or Cygwin) that don't provide this function, we implement our own version
   (e.g. mac/fmemopen.cpp).
 */
FILE *fmemopen (void *buf, size_t size, const char *opentype);

#endif

#endif
