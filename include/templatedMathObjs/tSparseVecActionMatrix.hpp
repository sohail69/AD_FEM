#pragma once
#include "../UtilityObjects/macros.hpp"
#include "../UtilityObjects/lowLevelMFEM.hpp"
#include "mfem.hpp"


//
// A sparse matrix implementation
// that takes is used purely for its
// functionality to multiply against
// a vector, this uses the CSR format
// it is made specifically for element
// interpolators
//
template<typename Numeric>
struct tSparseVecActionMat{
  //Variable vector and iterator
  mfem::Memory<int> offsets;
  mfem::Memory<int> columns;
  mfem::Memory<Numeric> Values;
  int size;


  //Matrix-Vector multiplication
  //for multi-vector support
  template<typename VecType1, typename VecType2>
  void Mult(const VecType1 x, VecType2 y);
};
