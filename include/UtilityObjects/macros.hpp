#pragma once
#include <array>
#include "lowLevelMFEM.hpp"

//The number definitions
using REAL64 = double;
using UINT64 = long unsigned;
using INT64  = long int;

//Inline the function
#define FORCE_INLINE MFEM_HOST_DEVICE inline __attribute__((always_inline))

//Pack the struct
#define PACKSTRUCT __attribute__ ((packed))
