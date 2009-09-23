#include "../compat.h"

#ifndef HAVE_FMEMOPEN

/* OS X and Cygwin don't have the fmemopen() function, so we provide our own, which just
 * saves the buffer to a file and then returns its file handle
 */


FILE *fmemopen (void *buf, size_t size, const char *opentype)
{
  FILE *f;
  assert(strcmp(opentype, "r") == 0);
  f = tmpfile();
  fwrite(buf, 1, size, f);
  rewind(f);
  return f;
}

#endif
