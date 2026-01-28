#pragma once
#include "../UtilityObjects/macros.hpp"
#include "../UtilityObjects/lowLevelMFEM.hpp"
#include "mfem.hpp"


/*****************************************\
!
!  This generates a single restriction 
!  operator that acts over all the Vars
!  to give the element vectors.
!
\*****************************************/
class tRestrictOperator : public mfem::Operator
{
  private:

  public:
    //Constructor
    tRestrictOperator();

    //Get the operator
    mfem::Operator & GetOperator() const;
};
