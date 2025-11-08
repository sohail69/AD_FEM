#pragma once
#include "../macros.hpp"

/***************************************\
!
!  Sparse Iterator data
!
\***************************************/
struct PACKSTRUCT iteratorData{
  FORCE_INLINE dualNumber(value_t r=0.0, gradient_t eps=0.0): val(r), grad(eps){};
};

