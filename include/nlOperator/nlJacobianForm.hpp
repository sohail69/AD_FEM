#pragma once
#include "mfem.hpp"
#include "../../../MFEM_STUFF/mfem-4.7/build/include/mfem/fem/intrules.hpp"
#include "../../../MFEM_STUFF/mfem-4.7/build/include/mfem/fem/lininteg.hpp"
#include "TCoeff.hpp"

using namespace std;
using namespace mfem;

/*****************************************\
!
! Jacobian integrator acts as a linearForm
!
!  Min[e(u_1,..,u_n)]
!  A.u = c
!
! The non-linear problem is solver by 
! first making the problem stationary with
! regards to the displacements:
! R_m(u) = de/du_m
!
/*****************************************/
template<typename Number>
class JacobianIntegrator : public LinearFormIntegrator
{
private:
  //Integration level

  //Energy Functional form coefficient
  TCoefficient<Number> *EFuncCoeff;

  //Reference to block vector of all data
  Array<int> btoffs;
  mutable BlockVector & xBlock;

  //Reference to block vector of all data
  Array<int> elm_btoffs;
  mutable BlockVector & elm_xBlock;
  DenseMatrix dshape, dshapedxt;
public:
   /// Constructs the domain integrator (Q, grad v)
   residualIntegrator(): LinearFormIntegrator(){}

   //Checks for device support
   bool SupportsDevice() const override { return true; }

   /** Given a particular Finite Element and a transformation (Tr)
       computes the element right hand side element vector, elvect. */
   void AssembleRHSElementVect(const FiniteElement &el,
                               ElementTransformation &Tr,
                               Vector &elvect) override;

   using LinearFormIntegrator::AssembleRHSElementVect;
};




/*****************************************\
!
! Residual integrator acts on a 
! RHSElement vector
!
/*****************************************/
template<typename Number>
void residualIntegrator<Number>::AssembleRHSElementVect(const FiniteElement &el
                                                      , ElementTransformation &Tr
                                                      , Vector &elvect)
{



};



