#pragma once
#include <functional>
#include "mfem.hpp"

/*****************************************\
!
!  Abstract input preprocessor used for
!  making input vector discrete data to
!  continuous function sampled data for
!  usage in the coefficients
!
\*****************************************/
template<typename Number>
class AbstractInputCleaner
{
  public:
   /// Define a time-independent templated coefficient
   AbstractInputCleaner(Array<int> blocks){};

 
};


/*****************************************\
!
!  Templated energy functional coefficient
!  integrator Base class:
!    Takes in a block vector and integrates
!    the energy functional.
!
\*****************************************/
template<typename Number>
class TCoefficientIntegrator
{
  public:
   /// Define a time-independent templated coefficient
   TCoefficientIntegrator(Array<int> used_blocks);

   /// Evaluate the integral of the element
   Number Eval(Array<int> elm_btoffs, mfem::Vector elm_x)=0;

   /// Coefficicient destructor
   ~TCoefficientIntegrator();
};
