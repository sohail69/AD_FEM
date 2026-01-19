#pragma once
#include "../UtilityObjects/macros.hpp"
#include "../UtilityObjects/lowLevelMFEM.hpp"
#include "mfem.hpp"


//
// tensor-multi-iterator Functor used to
// access flattened muli-dimensional vars
//
template<unsigned RANK>
struct PACKSTRUCT multiIterator{
  const unsigned rank;
  unsigned Sizes[RANK];

  //Accessor multi-iter flattening
  FORCE_INLINE unsigned operator()(const unsigned Iter[]){
    unsigned L, I=0;
    #pragma unroll
    for(unsigned J=0; J<rank; J++){
      L=Iter[J];
      #pragma unroll
      for(unsigned K=0; K<(rank-J); K++) L = L*Sizes[K];
      I=I+L;
    }
    return I;
  };
};

template<unsigned RANK>
struct PACKSTRUCT DevMultiIterator{
  unsigned Sizes[RANK];

  //Accessor multi-iter flattening
  FORCE_INLINE unsigned operator()(const unsigned Iter[]){
    unsigned L, I=0;
    #pragma unroll
    for(unsigned J=0; J<RANK; J++){
      L=Iter[J];
      #pragma unroll
      for(unsigned K=0; K<(RANK-J); K++) L = L*Sizes[K];
      I=I+L;
    }
    return I;
  };
};


//
// templated vector of numbers 
// using MFEMs memory layout
// with multiple iterator functors
//
template<typename Numeric>
struct tVarVectorMFEM{
  //Variable vector and iterator
  mfem::Memory<int> offsets;
  mfem::Memory<multiIterator<5>> Iters;
  mfem::Memory<Numeric> varData;
  int size;

  //Constructor and destructor
  tVarVectorMFEM(int size_, mfem::MemoryType mt): varData(size_, mt), size(size_){};

  ~tVarVectorMFEM(){varData.Delete();};

  //Data accessor
  FORCE_INLINE Numeric &operator()(unsigned VarID, unsigned Ip, unsigned Iter[]) {
    return varData[offsets[VarID] + Iters[VarID](Iter)];
  };
};








