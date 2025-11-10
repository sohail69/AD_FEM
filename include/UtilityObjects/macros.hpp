#pragma once
#include <omp.h>
#include <array>

//The number definitions
using REAL64 = double;
using UINT64 = long unsigned;
using INT64  = long int;

//Inline the function
#define FORCE_INLINE inline __attribute__((always_inline))

//Pack the struct
#define PACKSTRUCT __attribute__ ((packed))
