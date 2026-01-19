#pragma once
#include "../UtilityObjects/macros.hpp"
#include "../UtilityObjects/lowLevelMFEM.hpp"
#include "mfem.hpp"


/*****************************************\
!
!  This generates a single interpolator
!  sparse matrix from multiple fe-spaces
!  grid functions etc... to convert from
!  discrete data to sampled continuous data
!
\*****************************************/
template<typename Numeric>
class tInterpolator{


  void GenerateSparseOperator();
};
