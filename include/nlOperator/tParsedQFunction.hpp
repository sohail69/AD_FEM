#pragma once
#include <functional>
#include "../UtilityObjects/macros.hpp"
#include "../UtilityObjects/lowLevelMFEM.hpp"
#include "../templatedMathObjs/tMultiVarVector.hpp"
#include "mfem.hpp"

/*****************************************\
!
!  Templated energy functional coefficient
!  integrator Base class:
!   -Takes in a sampled Var and evaluates
!    the energy functional.
!
\*****************************************/
template<typename Number>
class tParsedQFunction
{
  private:
    //Integration rule ID of the coefficient
    //integrator. This class only takes in
    //functions and continuous variables,
    //it is not responsible for its own
    //integration
    unsigned IntegRuleID=0;

  public:
   /// Define a time-independent templated coefficient
   tParsedQFunction(Array<int> used_blocks, unsigned integID){};

   /// Coefficicient destructor
   ~tParsedQFunction(){};

   /// Evaluate the integral of the element
   unsigned &GetIntegRule(){return IntegRuleID;};

   /// Evaluate the integral of the element
   Number Eval(Array<int> InputBlocks, tVarVectorMFEM<Number> elm_vars){return 0;};
};
