/*
   This header file contains platform compatibility layer.

   It is included from common.h, so it is automatically included in all hermes
   sources. The implementation of the functions in this file is in the
   src/compat directory.
*/

#ifndef __HERMES2D_COMPAT_H
#define __HERMES2D_COMPAT_H


#ifndef HAVE_FMEMOPEN
/// Implementation of GNU fmemopen. Intended to be used if the current platform does not support it. 
FILE *fmemopen (void *buf, size_t size, const char *opentype);
#endif

//Windows DLL export/import definitions
#if defined(WIN32) || defined(_WINDOWS)
# if defined(_HERMESDLL)
#   define PUBLIC_API __declspec(dllexport)
#   define PUBLIC_API_USED_TEMPLATE(__implementation) template class PUBLIC_API __implementation
#   define PUBLIC_API_USED_STL_VECTOR(__type) PUBLIC_API_USED_TEMPLATE(std::allocator<__type>); PUBLIC_API_USED_TEMPLATE(std::vector<__type>)
#   define EXTERN extern PUBLIC_API
# else
#   define PUBLIC_API __declspec(dllimport)
#   define PUBLIC_API_USED_TEMPLATE(__implementation)
//#   define PUBLIC_API_USED_TEMPLATE(__implementation) extern template class PUBLIC_API __implementation
#   define PUBLIC_API_USED_STL_VECTOR(__type) PUBLIC_API_USED_TEMPLATE(std::allocator<__type>); PUBLIC_API_USED_TEMPLATE(std::vector<__type>)
#   define EXTERN extern PUBLIC_API
# endif
#else
# define PUBLIC_API
# define PUBLIC_API_USED_TEMPLATE(__implementation)
# define PUBLIC_API_USED_STL_VECTOR(__type)
# define EXTERN extern
#endif

//C99 functions
#include "compat/c99_functions.h"

#endif
