#ifndef __HERMES2D_C99_FUNCTIONS_H
#define __HERMES2D_C99_FUNCTIONS_H

#ifdef IMPLELENT_C99

/* Definitions of C99 specification. Used in a case of MSVC 2008 and 
 * below because MSVC follows C++ rather than C
 */

// Not-a-number constant.
EXTERN const double NAN;

// functions
EXTERN double exp2(double x); ///< exp 2
EXTERN double log2(double x); ///< log 2
EXTERN double cbrt(double x); ///< cubic root

#endif /* IMPLELENT_C99 */ 

#endif
