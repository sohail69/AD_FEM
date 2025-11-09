#pragma once
#include "../macros.hpp"
#include "tPowFuncs.hpp"


// Returns the absolute value of
// a function, doesn't work for
// complex numbers
template<typename Number>
FORCE_INLINE Number abs(Number z){
  return  sqrt<Number>(z*z);
};


// Returns the absolute value of
// a function, doesn't work for
// complex numbers
template<typename Number>
FORCE_INLINE Number fabs(Number z){
  return  sqrt<Number>(z*z);
};


// Returns the multiply add function
// for single number types
template<typename Number>
FORCE_INLINE Number fma(Number x, Number y, Number z){
  return x*y + z;
};

