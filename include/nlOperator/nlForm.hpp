#pragma once
#include "mfem.hpp"
//#include "../../../MFEM_STUFF/mfem-4.7/build/include/mfem/fem/lininteg.hpp"
#include "../../../MFEM_STUFF/mfem-4.7/build/include/mfem/fem/intrules.hpp"
#include "../../../MFEM_STUFF/mfem-4.7/build/include/mfem/linalg/dtensor.hpp"
#include "../templatedMathObjs/dualNumber.hpp"
#include "../templatedMathObjs/tVector.hpp"
#include "../UtilityObjects/utilityFuncs.hpp"
#include "TCoeffInteg.hpp"
#include <vector>

using namespace mfem;
template<typename Num> using dualSymNum = dualNumber<Num,Num>;
template<typename Num> using ddualSymNum = dualNumber<dualNumber<Num,Num>,dualNumber<Num,Num>>;

/*****************************************\
!
! Non-linear Form
!  solving this equation:
!  Min[e(u_1,..,u_n)]
!  A.u = c
!
! The non-linear problem is solved by 
! first making the functional e(u) 
! stationary with regards to the
! displacements:
! R_m(u) = de/du_m
!
! Commonly Newtons method is used for
! solving the problem and the Jacobian
! is needed:
! J_mn(u) = dR_m(u)/du_n
!
! This form does not take integrators as
! inputs, rather it uses a coefficient to
! directly derive the the residual and
! Jacobian forms
!
\*****************************************/
template<typename Number>
class nlForm : public Operator
{
private:
  //Reference to essential boundary conditions
  std::vector<ParFiniteElementSpace*> parFEs;
  std::vector<Array<int>*> ess_bcs_markers;
  Array<int> ess_bcs_tdofs;

  //Energy Functional form coefficient
  Array<TCoefficientIntegrator<Number>*> EFuncCoeff; //Coefficient Funcs
  std::vector<Array<int>*>               InpBlocks;  //Coefficient Input Vars

  //Total energy functionals
  std::function<dualSymNum<Number>(tVector<dualSymNum<Number>>)>  tEFunc1;
  std::function<ddualSymNum<Number>(tVector<ddualSymNum<Number>>)> tEFunc2;


  //Reference to block vector of element data
  Array<int> *btoffs_inp, *bvoffs_inp;
  Array<int> *btoffs_out, *bvoffs_out;
  mutable mfem::Vector  EBlockVector, EBlockResidual;
  mutable mfem::DenseMatrix elMat;

  //Used for directional derivatives Templated
  //dual number vector for Residual and Jacobian
  mutable tVector<dualSymNum<Number>>  *tExVec;
  mutable tVector<ddualSymNum<Number>> *tEJVec;
  mutable Operator *Jacobian_f=NULL;

  //Problem sizing and restriction operator
  const int nElms = 20, NDofs=8, NIntegs=1;
  const int n=nElms*NDofs;
  mfem::Operator *elem_restrict=NULL;

  //Device and memory configs
  const bool             & use_dev;
  const mfem::Device     & device;
  const mfem::MemoryType & mt;

public:
  //Constructor
  nlForm(const mfem::Device & device, const mfem::MemoryType & mt_, const bool & use_dev_);

  //Destructor
  ~nlForm();

  /// Adds an energy functional coeff to the
  /// list of them, these are differentiated
  /// to give the Residual form
  void AddEnergyFuncCoeff(TCoefficientIntegrator<Number> *EFuncCoeff, std::vector<int> InpBlocks);

  /// Assembles the residual form i.e. sums over all domain/bdr coefficients.
  /** When @ref UseFastAssembly "UseFastAssembly(true)" has been called and the
      assembly is compatible with device execution, it will be executed on the device. */
  virtual void Mult(const Vector & x, Vector & y) const;

  //Build the Jacobian Operator
  void buildJacobian(const Vector & x) const;

  //Returns a handle to the Jacobian
  mfem::Operator & GetGradient(const mfem::Vector &x) const override;
};


/*****************************************\
!
! Construct the residual form
!
\*****************************************/
template<typename Number>
nlForm<Number>::nlForm(const mfem::Device & dev, const mfem::MemoryType & mt_, const bool & usedev_):
                             Operator(n,n), device(dev), mt(mt_), use_dev(usedev_)
                           , EBlockVector(n,mt_), EBlockResidual(n,mt_)
{
  //////////////////////////
  ///Set sizes
  //////////////////////////
  elMat.SetSize(NDofs);
  EBlockVector.SetSize(n);
  EBlockResidual.SetSize(n);
  tExVec = new tVector<dualSymNum<Number>>(NDofs,mt);
  tEJVec = new tVector<ddualSymNum<Number>>(NDofs,mt);

  //////////////////////////
  ///Set the functions
  //////////////////////////
  tEFunc1;
  tEFunc2;


  //////////////////////////
  ///Check sizes
  //////////////////////////
  std::cout <<  elMat.Width()         << std::endl;
  std::cout <<  elMat.Height()        << std::endl;
  std::cout <<  tExVec->size          << std::endl;
  std::cout <<  tEJVec->size          << std::endl;
  std::cout <<  EBlockVector.Size()   << std::endl;
  std::cout <<  EBlockResidual.Size() << std::endl;
};


/*****************************************\
!
! Destroy the residual form
!
\*****************************************/
template<typename Number>
nlForm<Number>::~nlForm()
{
  delete tExVec, tEJVec;
};


/*****************************************\
!
!  (Build? and) return the Jacobian
!             operator
!
\*****************************************/
template<typename Number>
mfem::Operator & nlForm<Number>::GetGradient(const mfem::Vector &x) const
{
  buildJacobian(x);
  return *Jacobian_f;
};


/*****************************************\
!
!     Assemble the residual vector
!             and output
!
\*****************************************/
template<typename Number>
void nlForm<Number>::Mult(const Vector & x, Vector & y) const
{
  //Get the element data vectors
  if(elem_restrict != NULL) elem_restrict->Mult(x,EBlockVector);

  //Some stuff for devices
  const auto ElmVecs = Reshape(EBlockVector.Read(), NDofs, nElms);
//  auto ElmRess       = Reshape(EBlockResidual.ReadWrite(), NDofs, nElms);

  mfem::forall_switch(use_dev, nElms, [=] MFEM_HOST_DEVICE (int IElm)
  {
    //Copy element vector in
    for(int IDofs=0; IDofs<NDofs; IDofs++) (*tExVec)[IDofs].val = ElmVecs(IElm, IDofs);

    //Perturb each of the DOF's
    for(int IDofs=0; IDofs<NDofs; IDofs++){
      (*tExVec)[IDofs].grad = 1.0;
      EBlockResidual(IElm*NDofs + IDofs) += tEFunc1( *tExVec).grad;
      (*tExVec)[IDofs].grad = 0.0;
    }
  });

  //Get the residual vector
  if(elem_restrict        != NULL) elem_restrict->MultTranspose(EBlockResidual,y);
  if(ess_bcs_tdofs.Size() != 0)    y.SetSubVector(ess_bcs_tdofs,0.00);
};


/*****************************************\
!
!    Assemble and build the Jacobian
!               Operator
!
\*****************************************/
template<typename Number>
void nlForm<Number>::buildJacobian(const Vector & x) const
{
  //Get the element data vectors
  if(elem_restrict != NULL) elem_restrict->Mult(x,EBlockVector);
  const auto ElmVecs = Reshape(EBlockVector.Read(), NDofs, nElms);

  mfem::forall_switch(use_dev, nElms, [=] MFEM_HOST_DEVICE (int IElm)
  {
    //Copy element vector in
    for(int IDofs=0; IDofs<NDofs; IDofs++) (*tEJVec)[IDofs].val.val = ElmVecs(IDofs, IElm);

    //Perturb each of the DOF's in Ith-row
    for(int IDofs=0; IDofs<NDofs; IDofs++){
      (*tEJVec)[IDofs].val.grad = 1.0;
      //Perturb each of the DOF's in Jth-column
      for(int JDofs=0; JDofs<NDofs; JDofs++){
        (*tEJVec)[JDofs].grad.val = 1.0;
        elMat(IDofs, JDofs) += tEFunc2(*tEJVec).grad.grad;
        (*tEJVec)[JDofs].grad.val = 0.0;
      }
      (*tEJVec)[IDofs].val.grad = 0.0;
    }
  });
};
