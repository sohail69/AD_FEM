#include "../nlOperator/TCoeffInteg.hpp"


/*****************************************\
!
!  Example energy functional coefficients
!
\*****************************************/
template<typename Number>
class TLElasticityCoeffIntegrator
{
  public:
   /// Define a time-independent templated coefficient
   TCoefficientIntegrator(Array<int> used_blocks);

   /// Evaluate the integral of the element
   Number Eval(Array<int> elm_btoffs, BlockVector elm_x)=0;

   /// Coefficicient destructor
   ~TCoefficientIntegrator();
};

