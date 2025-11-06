#include "../nlOperator/TCoeffInteg.hpp"


/*****************************************\
!
!  Example energy functional coefficients
!
\*****************************************/
template<typename Number>
class TLElasticityCoeffIntegrator
{
  private:
    //consts// dShapeFunc, ShapeFunc, detJ, weights
    //Vars  // u, gradU 


/*
  //Element transforms etc..
  Array<int> vdofs;
  ElementTransformation *eltrans;
  DofTransformation *doftrans;
*/


  public:
   /// Define a time-independent templated coefficient
   TLElasticityCoeffIntegrator(Array<int> used_blocks);

   /// Evaluate the integral of the element
   Number Eval(Array<int> elm_btoffs, BlockVector elm_x)=0;

   /// Coefficicient destructor
   ~TLElasticityCoeffIntegrator();
};

