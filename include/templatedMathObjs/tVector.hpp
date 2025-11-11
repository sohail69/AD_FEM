#pragma once
#include "../UtilityObjects/macros.hpp"
#include "../UtilityObjects/lowLevelMFEM.hpp"
#include "mfem.hpp"

//
// templated vector of numbers 
//
template<typename Numeric>
struct tVector{
  Numeric *data;
  int size, Iter;

  tVector(){};
  tVector(const int size_, mfem::MemoryType mt): size(size_){data = new Numeric[size_];};
  ~tVector(){delete[] data;};

  //Access Data
  FORCE_INLINE Numeric &operator()(int I) {return data[I];};
  FORCE_INLINE Numeric &operator[](int I) {return data[I];};
};


//
// templated vector of numbers 
// using MFEMs memory layout
//
template<typename Numeric>
struct tVectorMFEM{
  //Definition of templated vector
  mfem::Memory<Numeric> data;
  int size, Iter;

  //Constructor and destructor
  tVectorMFEM(int size_, mfem::MemoryType mt): data(size_, mt), size(size_){};
  ~tVectorMFEM(){data.Delete();};

  //Access Data
  FORCE_INLINE Numeric &operator()(int I) {return data[I];};
  FORCE_INLINE Numeric &operator[](int I) {return data[I];};
};
