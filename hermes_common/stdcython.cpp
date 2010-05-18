#include <stdexcept>

#include "stdcython.h"

PyObject* global_empty_tuple;

void init_global_empty_tuple() {
  _CALLED_ONLY_ONCE;

  global_empty_tuple = PyTuple_New(0);
}

void throw_exception(char *text)
{
    throw std::runtime_error(text);
}
