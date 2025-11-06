#pragma once
#include "../macros.hpp"

// Returns the exponential
// of a number
template<typename Number>
FORCE_INLINE Number exp(const Number x)
{
  Number bigNum(33554432.0), one(1.0), xCopy(x);
  Number expX = one + (xCopy/bigNum);

  #pragma unroll
  for(unsigned I=0; I<25; I++) expX = expX*expX; // m=2^24
  return expX;
};


// Returns the natural log of a 
// number using 10 iterations
// of Halley's method
template<typename Number>
FORCE_INLINE Number log(const Number x)
{
  //B = log(x) => x = exp(B)
  Number B(0.00), expB(0.00), const2(2.0);
  Number xPb(0.00), xMb(0.00);

  //Use iterative formula x-exp(b) = 0
  //and solve with 5 iterations of
  //Halley's Method
  #pragma unroll
  for(unsigned I=0; I<10; I++){
    expB =  exp<Number>(B);
    xPb  = (expB + x);
    xMb  = (expB - x);
    B = B - const2*(xMb/xPb);
  }
  return B;
};


// Returns the common log
// of a number (base 10)
template<typename Number>
FORCE_INLINE Number log10(const Number x)
{
  Number ln10 = log<Number>(Number(10.0));
  Number lnX = log<Number>(x);
  return lnX/ln10;
};


// Returns the exponential function
// base 2
template<typename Number>
FORCE_INLINE Number exp2(const Number x)
{
  Number ln2 = log<Number>(Number(2.0));
  return exp<Number>(x*ln2);
};


// Returns the exponential function
// minus 1
template<typename Number>
FORCE_INLINE Number expm1(const Number x)
{
  Number one(1.0);
  return exp<Number>(x) - one;
};


// Returns the binary log
// of a number (base 2)
template<typename Number>
FORCE_INLINE Number log2(const Number x)
{
  Number ln2 = log<Number>(Number(2.0));
  Number lnX = log<Number>(x);
  return lnX/ln2;
};
