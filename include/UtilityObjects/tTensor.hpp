#pragma once
#include "macros.hpp"
#include <functional>


//Tensor Iterating function and
//inverse Iterating functions takes
//an array input
template<UINT64 DIM>
FORCE_INLINE UINT64 TensorIterator(const UINT64 Iters[DIM], const UINT64 Sizes[DIM])
{
 UINT64 I = 0;
  #pragma unroll
  for(UINT64 J=0; J<DIM; J++){
    UINT64 L = Iters[J];
    #pragma unroll
    for(UINT64 K=0; K<(DIM-J); K++){
      L = L*Sizes[J];
    }
    I = I + L;
  }
  return I;
};
/*
template<UINT64 DIM>
FORCE_INLINE void InvTensorIterator(const UINT64 & I, const UINT64 Sizes[DIM], UINT64 Iters[DIM])
{
  #pragma unroll
  for(UINT64 J=0; J<DIM; J++) Iters[J] = Sizes[J] - 1;
};*/

// Templated tensor 
template<typename Number, unsigned RANK>
struct PACKSTRUCT tTensor{
  Number *val=NULL;
  unsigned Size[RANK], tSize;

  //Constructor
  tTensor(const unsigned Size_[RANK])
  {
    for(unsigned I=0; I<RANK; I++){
     Size[I] = Size_[I];
     tSize  *= Size_[I];
    }
    if(tSize > 0) val = new Number[tSize];
  };

  //Destructor
  ~tTensor()
  {
    tSize = 0;
    if(val != NULL) delete[] val;
  };

  //Set-Retrieve data
  FORCE_INLINE Number & operator()(const UINT64 Iters[RANK]) 
  {
    unsigned K = TensorIterator<RANK>(Iters, Size);
    return val[K];
  };
};

//
// Matrix inversion special
// case tensor
//
template<typename Number>
void Invert2By2(const tTensor<Number,2> & Mat22, tTensor<Number,2> & invMat22)
{
  Number a, b, c, d, det;
  a = Mat22(0,0);  b = Mat22(0,1);
  c = Mat22(1,0);  d = Mat22(1,1);
  det = a*d - b*c;

  //Inversion
  invMat22(0,0) =  d/det;
  invMat22(0,1) = -b/det;
  invMat22(1,0) = -c/det;
  invMat22(1,1) =  a/det;
};


template<typename Number>
void Invert3by3(const tTensor<Number,2> & Mat33, tTensor<Number,2> & invMat33)
{
  Number a, b, c, d, e, f, g, h, i, det;
  a = Mat33(0,0);  b = Mat33(0,1);  c = Mat33(0,2);  
  d = Mat33(1,0);  e = Mat33(1,1);  f = Mat33(1,2);
  g = Mat33(2,0);  h = Mat33(2,1);  i = Mat33(2,2);
  det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);

  //Inversion
  invMat33(0,0) =  (e*i - f*h)/det;
  invMat33(0,1) = -(d*i - f*g)/det;
  invMat33(0,2) =  (d*h - e*g)/det;

  invMat33(1,0) = -(b*i - c*h)/det;
  invMat33(1,1) =  (a*i - c*g)/det;
  invMat33(1,2) = -(a*h - b*g)/det;

  invMat33(2,0) =  (b*f - c*e)/det;
  invMat33(2,1) = -(a*f - c*d)/det;
  invMat33(2,2) =  (a*e - b*d)/det;
};
