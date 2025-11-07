#pragma once
#include "mfem.hpp"

//The number definitions
using REAL64 = double;
using UINT64 = long unsigned;
using INT64  = long int;

//Inline the function
#define FORCE_INLINE MFEM_HOST_DEVICE inline __attribute__((always_inline))

//Pack the struct
#define PACKSTRUCT __attribute__ ((packed))

//parallel for loop
#define forAll(I, N) for(I=0; I<N; I++) 


//Tensor Iterating function and
//inverse Iterating functions takes
//an array input
template<UINT64 DIM>
FORCE_INLINE void TensorIterator(const UINT64 Iters[DIM], const UINT64 Sizes[DIM], UINT64 & I)
{
  I = 0;
  #pragma unroll
  for(UINT64 J=0; J<DIM; J++){
    UINT64 L = Iters[J];
    #pragma unroll
    for(UINT64 K=0; K<(DIM-J); K++){
      L = L*Sizes[J];
    }
    I = I + L;
  }
};

template<UINT64 DIM>
FORCE_INLINE void InvTensorIterator(const UINT64 & I, const UINT64 Sizes[DIM], UINT64 Iters[DIM])
{
  #pragma unroll
  for(UINT64 J=0; J<DIM; J++) Iters[J] = Sizes[J] - 1;
};
