#pragma once
#include "../macros.hpp"
#include "../../../../MFEM_STUFF/mfem-4.7/build/include/mfem/general/array.hpp"
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
  mfem::Memory<multiIterator> Iters;
  mfem::Memory<Numeric> varData;
  int size;

  //Constructor and destructor
  tVarVectorMFEM(int size_, mfem::MemoryType mt): data(size_, mt), size(size_){};
  ~tVarVectorMFEM(){data.Delete();};

  //Access Data
//  FORCE_INLINE Numeric &operator()(int I) {return data[I];};
//  FORCE_INLINE Numeric &operator[](int I) {return data[I];};
};


//
// multi-iterator used to
// access flattened muli-dimensional
// vars
//
struct multiIterator{
  unsigned size, Iter;

  //Constructor and destructor
  FORCE_INLINE multiIterator(const unsigned size, ){};

  //Accessor multi-iter flattening
  FORCE_INLINE unsigned operator()(const unsigned Iter[size]){return unsigned(0) };
};
