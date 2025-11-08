#pragma once
#include "../macros.hpp"

/***************************************\
!
!  Sparse Iterator data
!
\***************************************/
struct PACKSTRUCT iteratorData{
  int dim, iterStart, tSize;
  int *DimSizes;

//  FORCE_INLINE dualNumber(): val(r), grad(eps){};
};

