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
class tInterpolator : public mfem::DenseMatrix
{
  private:
    const mfem::Array<mfem::ParGridFunctions*> & ParGFS;
    mfem::DenseMatrix IOp;

  public:
    //Constructor
    tInterpolator();

    //Get the interpolation matrix
    mfem::DenseMatrix & GetMat(){return IOp};
};
