#pragma once
#include "../UtilityObjects/macros.hpp"
#include "../UtilityObjects/lowLevelMFEM.hpp"
#include "mfem.hpp"


/*****************************************\
!
!  This generates a single restriction 
!  operator that acts over all the Vars
!  to give the element vectors. (This is
!  for the total mesh level local 
!  restrictions on sub topologies on 
!  multi-meshes are handled internally by
!  the FE-spaces) 
!
\*****************************************/
template<typename UINT>
class tRestrictOperator : public mfem::Operator
{
  private:
    const mfem::Array<mfem::ParFiniteElementSpace*> & ParFEs;
    mfem::Array<UINT> row_offsets, column_offsets;
    mfem::BlockOperator *blockRestriction;

    UINT OperatorSizeM(const mfem::Array<mfem::ParFiniteElementSpace*> & ParFEs_);
    UINT OperatorSizeN(const mfem::Array<mfem::ParFiniteElementSpace*> & ParFEs_);
  public:
    //Constructor
    tRestrictOperator(const mfem::Array<mfem::ParFiniteElementSpace*> & ParFEs_);

    //Destructor
    ~tRestrictOperator();

    //Get the operator
    void Mult(const mfem::Vector & x, mfem::Vector & y) const;

    //Get the operator
    mfem::Operator & GetGradient() const;
};


/*****************************************\
!
!  This implements the tRestrictOperator
!  class
!
\*****************************************/
// The operator sizingfunction
// for MxN dimensions
template<typename UINT>
UINT tRestrictOperator<UINT>::OperatorSizeM(const mfem::Array<mfem::ParFiniteElementSpace*> & ParFEs_)
{
  UINT EDofs=0;
  for(UINT I=0; I<ParFEs_->Size(); I++ ) EDofs += ParFEs_[I]->GetVSize();
  return EDofs;
};

template<typename UINT>
UINT tRestrictOperator<UINT>::OperatorSizeN(const mfem::Array<mfem::ParFiniteElementSpace*> & ParFEs_)
{
  UINT TDofs=0; 
  for(UINT I=0; I<ParFEs_->Size(); I++ ) TDofs += ParFEs_[I]->GlobalTrueVSize();
  return TDofs;
};

// The constructor
template<typename UINT>
tRestrictOperator<UINT>::tRestrictOperator(const mfem::Array<mfem::ParFiniteElementSpace*> & ParFEs_):
                                           mfem::Operator(OperatorSizeM(ParFEs_),OperatorSizeN(ParFEs_))
                                         , row_offsets(ParFEs_.Size()), column_offsets(ParFEs_.Size())
                                         , ParFEs(ParFEs_)
{
  for(UINT I=0; I<ParFEs_->Size(); I++ ) row_offsets[I] = ParFEs_[I]->GetVSize();
  for(UINT I=0; I<ParFEs_->Size(); I++ ) column_offsets[I] = ParFEs_[I]->GlobalTrueVSize();
  blockRestriction = new mfem::BlockOperator(row_offsets, column_offsets);

  for(UINT I=0; I<ParFEs_->Size(); I++){
    //For diagonal problems
//    blockRestriction->SetBlock(I,I, ParFEs_.Ptr());
    //For off-diagonal problems
    for(UINT J=0; J<ParFEs_->Size(); J++){
      if(I != J){
//        if( ) blockRestriction->SetBlock(I,J, ParFEs_.Ptr());
      }
    }
  }
};


// The destructor
template<typename UINT>
tRestrictOperator<UINT>::~tRestrictOperator()
{
  delete blockRestriction;
};

// Mult (Apply the restriction operator)
template<typename UINT>
void tRestrictOperator<UINT>::Mult(const mfem::Vector & x, mfem::Vector & y) const
{
 blockRestriction->Mult(x,y);
};

// GetGradient (get the restriction Operator)
template<typename UINT>
mfem::Operator & tRestrictOperator<UINT>::GetGradient() const
{
  return *blockRestriction;
};
