// Triangle does not give C++-compatibility macros, so include this header instead.
// This also includes the triangle header in library-mode, and specifies the REAL and VOID types.
#ifndef TRIANGLE_LIBRARY_HEADER_WRAPPER_H

typedef double REAL;
typedef void VOID;

#ifdef __cplusplus
extern "C" {
#endif

#ifndef ANSI_DECLARATORS
#define ANSI_DECLARATORS
#endif
#ifndef TRILIBRARY
#define TRILIBRARY
#endif
#include "triangle/triangle.h"

#ifdef __cplusplus
} // end extern "C"
#endif

#define TRIANGLE_LIBRARY_HEADER_WRAPPER_H
#endif // TRIANGLE_LIBRARY_HEADER_WRAPPER_H
