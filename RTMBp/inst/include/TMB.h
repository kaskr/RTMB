#define TMB_SKINNY
// TMB already included => skip!
#ifndef TMB_OBJECTIVE_PTR

#include "config.h"
#ifndef TMB_H
#define TMB_H
#ifdef TMB_PRECOMPILE
/** \file
    \brief Include this file to extract declarations, definitions and selected code for pre-compilation
*/
#undef WITH_LIBTMB
#undef TMB_PRECOMPILE
#undef CSKIP
#undef IF_TMB_PRECOMPILE
#undef TMB_EXTERN
// Redefine
#undef  WITH_LIBTMB
#define TMB_PRECOMPILE
#define CSKIP(...) __VA_ARGS__
#define IF_TMB_PRECOMPILE(...) __VA_ARGS__
#define TMB_EXTERN
#else
/** \file
    \brief Include this file to extract declarations only
*/
#undef WITH_LIBTMB
#undef TMB_PRECOMPILE
#undef CSKIP
#undef IF_TMB_PRECOMPILE
#undef TMB_EXTERN
// Redefine
#define WITH_LIBTMB
#undef  TMB_PRECOMPILE
#define CSKIP(...) ;
#define IF_TMB_PRECOMPILE(...)
#define TMB_EXTERN extern
#endif
#include <TMB.hpp>
#endif

#endif
