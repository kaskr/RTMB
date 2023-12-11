// Set the 'skinny' flag - otherwise same as RTMB/src/TMB.h
#define TMB_SKINNY
#include "config.h"
#ifndef TMB_H
#define TMB_H
#ifndef TMB_PRECOMPILE
#define WITH_LIBTMB
#endif
#include <TMB.hpp>
#endif

void rtmb_set_shared_pointers();
