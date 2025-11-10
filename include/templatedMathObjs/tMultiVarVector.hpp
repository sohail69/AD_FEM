#pragma once
#include "../macros.hpp"
#include "../../../../MFEM_STUFF/mfem-4.7/build/include/mfem/general/globals.hpp"
#include "../../../../MFEM_STUFF/mfem-4.7/build/include/mfem/general/mem_manager.hpp"
#include "../../../../MFEM_STUFF/mfem-4.7/build/include/mfem/general/device.hpp"

#include "mfem.hpp"

//
// templated vector of numbers 
// using MFEMs memory layout
// with multiple iterator functors
//

template<typename Numeric>
struct tVarVectorMFEM{
  //Definition of templated vector
  mfem::Memory<int> offsets;
//  mfem::Memory<multiIterator> Iters;
  mfem::Memory<Numeric> varData;
  int size;

  //Constructor and destructor
  tVarVectorMFEM(int size_, mfem::MemoryType mt): varData(size_, mt), size(size_){};
  ~tVarVectorMFEM(){varData.Delete();};

  //Access Data
//  FORCE_INLINE Numeric &operator()(int I) {return varData[I];};
//  FORCE_INLINE Numeric &operator[](int I) {return varData[I];};
};


//
// multi-iterator used to
// access flattened muli-dimensional
// vars
//
template<unsigned Dim>
struct multiIterator{
  static const unsigned dim=Dim;
  unsigned *Sizes, *Iter;

  //Constructor and destructor
  multiIterator(){ Iter=new unsigned[Dim]; Sizes=new unsigned[Dim];};
  ~multiIterator(){ delete[] Iter, Sizes;};

  //Accessor multi-iter flattening
  FORCE_INLINE unsigned operator()(const unsigned Iter[dim]){return 0; };
};
