#pragma once
#include "../macros.hpp"

//
// Dual numbers 
//
template<typename real_t, typename imaginary_t>
struct PACKSTRUCT complexNumber{
  //Definition of dual number
  real_t      real;
  imaginary_t imag;
  complexNumber(real_t r = 0.0, imaginary_t i = 0.0) : real(r), imag(i) {}

  /***************************************\
  !
  !  Complex-Complex number operations
  !
  \***************************************/
  //Multiplication operator
  FORCE_INLINE complexNumber operator*(const complexNumber other) 
  {
    complexNumber newNum;
    newNum.real = real*other.real - imag*other.imag;
    newNum.imag = real*other.imag + imag*other.real;
    return newNum;
  };

  //Division operator
  FORCE_INLINE complexNumber operator/(const complexNumber other) 
  {
    complexNumber newNum;
    newNum.real = (real*other.real + imag*other.imag)/(other.real*other.real + other.imag*other.imag);
    newNum.imag = (imag*other.real - real*other.imag)/(other.real*other.real + other.imag*other.imag);
    return newNum;
  };

  //Addition operator
  FORCE_INLINE complexNumber operator+(const complexNumber other)
  {
    complexNumber newNum;
    newNum.real = real + other.real;
    newNum.imag = imag + other.imag;
    return newNum;
  };

  //Subtraction operator
  FORCE_INLINE complexNumber operator-(const complexNumber other) 
  {
    complexNumber newNum;
    newNum.real = real - other.real;
    newNum.imag = imag - other.imag;
    return newNum;
  };
};


/***************************************\
!
!  Number-Complex number operations
!
\***************************************/
//Multiplication operator
template<typename val_t, typename imag_t>
FORCE_INLINE complexNumber<val_t,imag_t> operator*(const complexNumber<val_t,imag_t> dNum, const double Num) 
{
  complexNumber<val_t,imag_t> newVal;
  newVal.real = Num*dNum.real;
  newVal.imag = Num*dNum.imag;
  return newVal;
};

template<typename val_t, typename imag_t>
FORCE_INLINE complexNumber<val_t,imag_t> operator*(const complexNumber<val_t,imag_t> dNum, const float Num) 
{
  complexNumber<val_t,imag_t> newVal;
  newVal.real = Num*dNum.real;
  newVal.imag = Num*dNum.imag;
  return newVal;
};

template<typename val_t, typename imag_t, typename Number>
FORCE_INLINE complexNumber<val_t,imag_t> operator*(Number Num, complexNumber<val_t,imag_t> dNum) 
{
  return dNum*Num;
};



