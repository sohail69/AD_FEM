#pragma once
#include "mfem.hpp"
//#include "../../../MFEM_STUFF/mfem-4.7/build/include/mfem/fem/lininteg.hpp"
#include "../../../MFEM_STUFF/mfem-4.7/build/include/mfem/fem/intrules.hpp"
#include "../../../MFEM_STUFF/mfem-4.7/build/include/mfem/linalg/dtensor.hpp"
#include "../templatedMathObjs/dualNumber.hpp"
#include "../templatedMathObjs/tVector.hpp"
#include "../UtilityObjects/utilityFuncs.hpp"
#include "TCoeffInteg.hpp"

using namespace mfem;

/*****************************************\
!
! Residual Form
!  solving this equation:
!  Min[e(u_1,..,u_n)]
!  A.u = c
!
! The non-linear problem is solved by 
! first making the problem stationary with
! regards to the displacements:
! R_m(u) = de/du_m
!
! This form does not take integrators as
! inputs, rather it uses a coefficient to
! directly derive the the residual
!
/*****************************************/
template<typename Number>
class residualForm : public Vector
{
private:
  //Energy Functional form coefficient
  Array<TCoefficientIntegrator<Number>*>       EFuncCoeff; //Coefficient Funcs
  std::vector<Array<int>*>           InpBlocks;  //Coefficient Input Vars
  tVector<dualNumber<Number,Number>> ElmDualData;//Dual vector containing element data

  //Reference to block vector of element data
  Array<int> *btoffs_inp, *bvoffs_inp;
  Array<int> *btoffs_out, *bvoffs_out;
//  mutable BlockVector xBlockElm, ;
  mfem::Vector                       EBlockVector, EBlockResidual;
  tVector<dualNumber<Number,Number>> tEVec;        //Templated dual number vector for auto-diff



  //The device and sizing based things
  bool use_dev=true;
  const UINT64 nElms = 20, NDofs=8, NIntegs=0;
  const UINT64 n = nElms*NDofs;
  const char *device_config = "cpu";
  const mfem::Device     & device;
  const mfem::MemoryType & mt;

  //Preprocesses element DOF-data
  //into continuous data used for coefficients
  void PreprocessData(const Vector & elm_x){};

public:
  //Constructor
  residualForm(const mfem::Device & device, const mfem::MemoryType & mt_);

   /// Adds an energy functional coeff to the
   /// list of them, these are differentiated
   /// to give the Residual form
   void AddEnergyFuncCoeff(TCoefficientIntegrator<Number> *EFuncCoeff, std::vector<int> InpBlocks);

   /// Assembles the residual form i.e. sums over all domain/bdr coefficients.
   /** When @ref UseFastAssembly "UseFastAssembly(true)" has been called and the
       assembly is compatible with device execution, it will be executed on the device. */
   void Assemble(const Vector & x);

   /// Assemble the vector on the true dofs, i.e. P^t v.
   void ParallelAssemble(Vector & x);
};

/*****************************************\
!
! Construct the residual form
!
/*****************************************/
template<typename Number>
residualForm<Number>::residualForm(const mfem::Device & device_, const mfem::MemoryType & mt_) :
 device(device_), mt(mt_), tEVec(NDofs,mt_), EBlockVector(n,mt_), EBlockResidual(n,mt_)
{};

/*****************************************\
!
! Assemble the residual vector
! in an element by element approach
!
/*****************************************/
template<typename Number>
void residualForm<Number>::Assemble(const Vector & x)
{
  //Get the element block vectors
//  gather(x,xBlockElm);


  //Some stuff for devices
  auto ElmVecs = Reshape(EBlockVector.Read(use_dev), nElms, NDofs);
  auto ElmRess = Reshape(EBlockResidual.ReadWrite(use_dev),  nElms, NDofs);
  auto d_X = tEVec.ReadWrite(use_dev);

  mfem::forall_switch(use_dev, nElms, [=] MFEM_HOST_DEVICE (int IElm)
  {
    //Copy element vector in
    for(UINT64 IDofs=0; IDofs<NDofs; IDofs++) d_X[IDofs].val = ElmVecs(IElm,IDofs);

    //Perturb each of the DOF's
    for(UINT64 IDofs=0; IDofs<NDofs; IDofs++){
      d_X[IDofs].grad = 1.0;
      //Add on each of the integrator contribution
      for(UINT64 IIntegs=0; IIntegs<NIntegs; IIntegs++){
        ElmRess(IElm,IDofs) +=  d_X[IDofs].grad;
      }
      d_X[IDofs].grad = 0.0;
    }

    //Copy element vector out
//    for(UINT64 IDofs=0; IDofs<NDofs; IDofs++) ElmRess(IElm,IDofs) = d_X[IDofs].grad;
  });
};
