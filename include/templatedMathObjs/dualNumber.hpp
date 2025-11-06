#pragma once
#include "../macros.hpp"

/***************************************\
!
!  Dual-Numbers
!
\***************************************/
template<typename value_t, typename gradient_t>
struct PACKSTRUCT dualNumber{
  //Definition of dual number
  value_t     val;
  gradient_t  grad;

  dualNumber(value_t r=0.0, gradient_t eps=0.0) : val(r), grad(eps){};
};

/***************************************\
!
!  Dual-Dual number operations
!
\***************************************/
//Multiplication operator
template<typename v_t, typename g_t>
FORCE_INLINE dualNumber<v_t,g_t> operator*(const dualNumber<v_t,g_t> a, const dualNumber<v_t,g_t> b) 
{
  dualNumber<v_t,g_t> newVal;
  newVal.val  = a.val*b.val;
  newVal.grad = a.val*b.grad + a.grad*b.val;
  return newVal;
};

//Division operator
template<typename v_t, typename g_t>
FORCE_INLINE dualNumber<v_t,g_t> operator/(const dualNumber<v_t,g_t> a, const dualNumber<v_t,g_t> b) 
{
  dualNumber<v_t, g_t> newVal;
  newVal.val  = (a.val/b.val);
  newVal.grad = (a.grad*b.val + a.val*b.grad)/(b.val*b.val);
  return newVal;
};

//Addition operator
template<typename v_t, typename g_t>
FORCE_INLINE dualNumber<v_t,g_t> operator+(const dualNumber<v_t,g_t> a, const dualNumber<v_t,g_t> b)
{
  dualNumber<v_t, g_t> newVal;
  newVal.val  = a.val + b.val;
  newVal.grad = a.grad + b.grad;
  return newVal;
};

//Subtraction operator
template<typename v_t, typename g_t>
FORCE_INLINE dualNumber<v_t,g_t> operator-(const dualNumber<v_t,g_t> a, const dualNumber<v_t,g_t> b)
{
  dualNumber<v_t, g_t> newVal;
  newVal.val  = a.val - b.val;
  newVal.grad = a.grad - b.grad;
  return newVal;
};


/***************************************\
!
!  Number-Dual number operations
!
\***************************************/
///////////
//Multiplication operator
///////////
template<typename val_t, typename grad_t>
FORCE_INLINE dualNumber<val_t,grad_t> operator*(const dualNumber<val_t,grad_t> dNum, const double Num) 
{
  dualNumber<val_t,grad_t> newVal;
  newVal.val  = Num*dNum.val;
  newVal.grad = Num*dNum.grad;
  return newVal;
};

template<typename val_t, typename grad_t>
FORCE_INLINE dualNumber<val_t,grad_t> operator*(const dualNumber<val_t,grad_t> dNum, const float Num) 
{
  dualNumber<val_t,grad_t> newVal;
  newVal.val  = Num*dNum.val;
  newVal.grad = Num*dNum.grad;
  return newVal;
};

template<typename val_t, typename grad_t, typename Number>
FORCE_INLINE dualNumber<val_t,grad_t> operator*(const Number Num, const dualNumber<val_t,grad_t> dNum) 
{
  return dNum*Num;
};

///////////
//Division operator
///////////
template<typename val_t, typename grad_t>
FORCE_INLINE dualNumber<val_t,grad_t> operator/(const dualNumber<val_t,grad_t> dNum, const double Num)
{
  dualNumber<val_t,grad_t> newVal;
  newVal.val  = newVal.val/Num;
  newVal.grad = newVal.grad/Num;
  return newVal;
};

template<typename val_t, typename grad_t>
FORCE_INLINE dualNumber<val_t,grad_t> operator/(const dualNumber<val_t,grad_t> dNum, const float Num)
{
   dualNumber<val_t,grad_t> newVal;
  newVal.val  = newVal.val/Num;
  newVal.grad = newVal.grad/Num;
  return newVal;
};

template<typename val_t, typename grad_t, typename Number>
FORCE_INLINE dualNumber<val_t,grad_t> operator/(const Number Num, const dualNumber<val_t,grad_t> dNum)
{
  dualNumber<val_t,grad_t> newVal(Num);
  return Num/dNum;
};

///////////
//Addition operator
///////////
template<typename val_t, typename grad_t>
FORCE_INLINE dualNumber<val_t,grad_t> operator+(const dualNumber<val_t,grad_t> dNum, const double Num)
{
  dualNumber<val_t,grad_t> newVal;
  newVal.val  = dNum.val  + Num;
  newVal.grad = dNum.grad;
  return newVal;
};

template<typename val_t, typename grad_t>
FORCE_INLINE dualNumber<val_t,grad_t> operator+(const dualNumber<val_t,grad_t> dNum, const float Num)
{
  dualNumber<val_t,grad_t> newVal;
  newVal.val  = dNum.val  + Num;
  newVal.grad = dNum.grad;
  return newVal;
};

template<typename val_t, typename grad_t, typename Number>
FORCE_INLINE dualNumber<val_t,grad_t> operator+(const Number Num, const dualNumber<val_t,grad_t> dNum)
{
  dualNumber<val_t,grad_t> newVal(Num);
  return dNum + Num;
};

///////////
//Subtraction operator
///////////
template<typename val_t, typename grad_t>
FORCE_INLINE dualNumber<val_t,grad_t> operator-(const dualNumber<val_t,grad_t> dNum, const double Num)
{
  dualNumber<val_t,grad_t> newVal;
  newVal.val  = dNum.val  - Num;
  newVal.grad = dNum.grad;
  return newVal;
};

template<typename val_t, typename grad_t>
FORCE_INLINE dualNumber<val_t,grad_t> operator-(const dualNumber<val_t,grad_t> dNum, const float Num)
{
  dualNumber<val_t,grad_t> newVal;
  newVal.val  = dNum.val  - Num;
  newVal.grad = dNum.grad;
  return newVal;
};

template<typename val_t, typename grad_t, typename Number>
FORCE_INLINE dualNumber<val_t,grad_t> operator-(const Number Num, const dualNumber<val_t,grad_t> dNum)
{
  dualNumber<val_t,grad_t> newVal(Num);
  return dNum - Num;
};
