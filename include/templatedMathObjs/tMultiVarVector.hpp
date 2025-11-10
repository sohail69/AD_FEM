#pragma once
#include "../UtilityObjects/macros.hpp"
#include "../../../../MFEM_STUFF/mfem-4.7/build/include/mfem/general/globals.hpp"
#include "../../../../MFEM_STUFF/mfem-4.7/build/include/mfem/general/mem_manager.hpp"
#include "../../../../MFEM_STUFF/mfem-4.7/build/include/mfem/general/device.hpp"
#include "mfem.hpp"


//
// tensor-multi-iterator Functor used to
// access flattened muli-dimensional vars
//
struct multiIterator{
  const unsigned rank;
//  unsigned Sizes[RANK];

  //Constructor
//  template<unsigned RANK>
//  multiIterator(unsigned Sizes_[RANK]):rank(RANK){
 //   #pragma unroll
//    for(unsigned J=0; J<rank; J++) Sizes[J]=Sizes_[J];
 // };

  //Accessor multi-iter flattening
  FORCE_INLINE unsigned operator()(const unsigned Iter[]){
    unsigned L, I=0;
    #pragma unroll
    for(unsigned J=0; J<rank; J++){
      L=Iter[J];
      #pragma unroll
//      for(unsigned K=0; K<(rank-J); K++) L = L*Sizes[K];
      I=I+L;
    }
//    return I;
    return 0;
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
  mfem::Memory<multiIterator> Iters;
  mfem::Memory<Numeric> varData;
  int size;

  //Constructor and destructor
  tVarVectorMFEM(int size_, mfem::MemoryType mt): varData(size_, mt), size(size_){};
  ~tVarVectorMFEM(){varData.Delete();};

  //Data accessor
  FORCE_INLINE Numeric &operator()(unsigned VarID, unsigned Iter[]) {
    return varData[offsets[VarID] + Iters[VarID](Iter)];
  };
};








