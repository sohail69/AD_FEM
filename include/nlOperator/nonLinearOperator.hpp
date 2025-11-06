#pragma once
#include "mfem.hpp"
#include "TCoeffInteg.hpp"

using namespace std;
using namespace mfem;

/*****************************************\
!
! Solves a general minimisation problem of
! non-linear energy functional:
!
!  Min[e(u_1,..,u_n)]
!  A.u = c
!
! The non-linear problem is solver by 
! first making the problem stationary with
! regards to the displacements:
! R_m(u) = de/du_m
!
! Commonly Newtons method is used for
! solving the problem and the Jacobian
! is needed:
! J_mn(u) = dR_m(u)/du_n
!
/*****************************************/
class nonLinearOperator : public Operator
{
protected:
   // FE-spaces and forms
   Array<ParFiniteElementSpace*> & fes;
   ParLinearForm   *Residual_f;
   ParBilinearForm *Jacobian_f;

   // Jacobian matrix storage
   mutable OperatorPtr opJac;

   //Coefficients
   ConstantCoefficient *lmbdaInv;

   //Internal gridfunctions used for coefficient
   //update and boundary conditions
   Array<int> btoffs;
   mutable Array<ParGridFunction*> gFuncs;
   mutable BlockVector xBlock, rhsBlock;

   //Dirchelet boundary condition
   Array<int> ess_bcs_markers;

   //Build the operators for the Jacobian
   void reassembleJacobian() const;
public:
   //Constructs fibre map operator
   nonLinearOperator(Array<ParFiniteElementSpace*> & fes_);

   //Calculates and returns the Residual
   virtual void Mult(const Vector &x, Vector &Residual) const;

   //Returns a handle to the Jacobian
   mfem::Operator & GetGradient(const mfem::Vector &x) const override;

   //Destroys the fibre map operator
   ~nonLinearOperator();
};


/*****************************************\
!
! Construct the problem Operator
!
\*****************************************/
nonLinearOperator::nonLinearOperator(Array<ParFiniteElementSpace*> & fes_): fes(fes_)
{


};


/*****************************************\
!
  //Build the operators for the Jacobian
!
\*****************************************/
void nonLinearOperator::reassembleJacobian() const
{

};


/*****************************************\
!
   //Calculates and returns the Residual
   //and recalculates the Jacobian
!
\*****************************************/
void nonLinearOperator::Mult(const Vector &x, Vector &Residual) const
{


};


/*****************************************\
!
   //Calculates and returns the Residual
   //and recalculates the Jacobian
!
\*****************************************/
 mfem::Operator & nonLinearOperator::GetGradient(const mfem::Vector &x) const
{

}


/*****************************************\
!
   //Destroys the fibre map operator
!
\*****************************************/
nonLinearOperator::~nonLinearOperator()
{

};




