#pragma once
#include "../macros.hpp"
#include "tExpLogFuncs.hpp"
#include "tPowFuncs.hpp"

// Returns the hyperbolic Cosine
// function of a number
template<typename Number>
Number cosh(Number theta){
  Number half(0.5), expX(exp<Number>(theta)), expMX(exp<Number>(-theta));
  return half*(expX + expMX);
};


// Returns the hyperbolic Sine
// function of a number
template<typename Number>
Number sinh(Number theta){
  Number half(0.5), expX(exp<Number>(theta)), expMX(exp<Number>(-theta));
  return half*(expX - expMX);
};


// Returns the hyperbolic tangent
// function of a number
template<typename Number>
Number tanh(Number theta){
  Number half(0.5), expX(exp<Number>(theta)), expMX(exp<Number>(-theta));
  return (expX - expMX)/(expX + expMX);
};


// Returns the area hyperbolic cosine
// function of a number
template<typename Number>
Number acosh(Number x){
  Number one(1.00);
  return  log<Number>(x + sqrt<Number>(x*x - one) );
};


// Returns the area hyperbolic sine
// function of a number
template<typename Number>
Number asinh(Number x){
  Number one(1.00);
  return log<Number>(x + sqrt<Number>(x*x + one) );
};


// Returns the area hyperbolic tangent
// function of a number
template<typename Number>
Number atanh(Number x){
  Number half(0.50), one(1.00);
  return half*log<Number>( (one + x)/(one - x) );
};
