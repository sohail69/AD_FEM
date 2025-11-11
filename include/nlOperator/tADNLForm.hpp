#pragma once
#include "mfem.hpp"
#include "../UtilityObjects/macros.hpp"
#include "../UtilityObjects/lowLevelMFEM.hpp"
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
class tADNLForm : public Operator
{
private:
  //Reference to essential boundary conditions
  const std::vector<ParGridFunction*> & Vars;
  std::vector<Array<int>*>      ess_bcs_markers;
  Array<int>                    ess_bcs_tdofs;

  //Energy Functional forms coefficient
  std::vector<TCoefficientIntegrator<Number>*>      EFuncCoeff; //Coefficient Funcs
  std::vector<std::function<dualSymNum<Number>()>>  dEfuncs;    //dual functions
  std::vector<std::function<ddualSymNum<Number>()>> ddEfuncs;   //dual dual functions

  //Reference to block vector of element data
  Array<int> *btoffs_inp, *bvoffs_inp;
  Array<int> *btoffs_out, *bvoffs_out;
  mutable mfem::Vector  *EBlockVector, *EBlockResidual;
  mutable mfem::DenseMatrix elMats;

  //Used for directional derivatives Templated
  //dual number vector for Residual and Jacobian
  mutable tVector<dualSymNum<Number>>  *tExVec;
  mutable tVector<ddualSymNum<Number>> *tEJVec;
  mutable Operator *Jacobian_f=NULL;

  //Used for Isogeometric interpolator functions
  //at all integration points for a given integration
  //rule
//  Array<IntegrationRule> IntRules;


  //Problem sizing and restriction operator
  const int NIntegs = 1;
  int nElms=0, NDofsMax=0, nElmDofs=0, NEQs=0;
  mfem::Operator *elem_restrict=NULL;

  //Device and memory configs
  const bool             & use_dev;
  const mfem::Device     & device;
  const mfem::MemoryType & mt;

public:
  //Constructor
  tADNLForm(const std::vector<ParGridFunction*> & Vars_, const mfem::Device & dev
          , const mfem::MemoryType & mt_, const bool & use_dev_);

  //Destructor
  ~tADNLForm();

  /// Adds an energy functional coeff to the
  /// list of them, these are differentiated
  /// to give the Residual form
  void AddEnergyFuncCoeff(TCoefficientIntegrator<Number> *EFuncCoeff, std::vector<int> InpBlocks);

  /// Assembles the residual form i.e.
  /// sums over all domain/bdr coefficients.
  virtual void Mult(const Vector & x, Vector & y) const;

  //Build the Jacobian Operator
  void buildJacobian(const Vector & x) const;

  //Returns a handle to the Jacobian
  mfem::Operator & GetGradient(const mfem::Vector &x) const override;
};


/*****************************************\
!
! Construct the non-linear form
!
\*****************************************/
template<typename Number>
tADNLForm<Number>::tADNLForm(const std::vector<ParGridFunction*> & Vars_, const mfem::Device & dev
                           , const mfem::MemoryType & mt_, const bool & use_dev_):
                             Vars(Vars_), device(dev), mt(mt_), use_dev(use_dev_)
{
  //////////////////////////
  ///Recover the object sizes
  //////////////////////////
  nElms = Vars[0]->ParFESpace()->GetMesh()->GetNE();
  for(int I=0; I<Vars.size(); I++){
    NEQs += Vars[I]->ParFESpace()->GetTrueVSize();
    int NDofVar=0;
    for(int J=0; J<nElms; J++){
      DofTransformation doftrans;
      Array<int> dofs;
      Vars[I]->ParFESpace()->GetElementDofs(J,dofs,doftrans);
      NDofVar = max(dofs.Size(),NDofVar);
      nElmDofs += dofs.Size();
    }
    NDofsMax += NDofVar;
  }

  //////////////////////////
  ///Allocate the memory objects
  //////////////////////////
  elMats.SetSize(NDofsMax);
  EBlockVector   = new mfem::Vector(nElmDofs,mt_);
  EBlockResidual = new mfem::Vector(nElmDofs,mt_);
  tExVec = new tVector<dualSymNum<Number>>(NDofsMax,mt);
  tEJVec = new tVector<ddualSymNum<Number>>(NDofsMax,mt);

  //////////////////////////
  ///Check sizes
  //////////////////////////
  std::cout <<  tExVec->size           << std::endl;
  std::cout <<  tEJVec->size           << std::endl;
  std::cout <<  EBlockVector->Size()   << std::endl;
  std::cout <<  EBlockResidual->Size() << std::endl;
};


/*****************************************\
!
! Destroy the Form
!
\*****************************************/
template<typename Number>
tADNLForm<Number>::~tADNLForm()
{
  delete tExVec, tEJVec, EBlockVector, EBlockResidual;
};


/*****************************************\
!
! Add energy functional coeffs
!
\*****************************************/
template<typename Number>
void tADNLForm<Number>::AddEnergyFuncCoeff(TCoefficientIntegrator<Number> *EFuncCoeff
                                         , std::vector<int> InpBlocks)
{


};


/*****************************************\
!
!  (Build? and) return the Jacobian
!             operator
!
\*****************************************/
template<typename Number>
mfem::Operator & tADNLForm<Number>::GetGradient(const mfem::Vector &x) const
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
void tADNLForm<Number>::Mult(const Vector & x, Vector & y) const
{
  //Get the element data vectors
  if(elem_restrict != NULL) elem_restrict->Mult(x,*EBlockVector);

  //Some stuff for devices
  const auto ElmVecs = Reshape(EBlockVector->Read(), NDofsMax, nElms);
  auto ElmRess       = Reshape(EBlockResidual->ReadWrite(), NDofsMax, nElms);

  mfem::forall_switch(use_dev, nElms, [=] MFEM_HOST_DEVICE (int IElm)
  {
    //Copy element vector in
    for(int IDofs=0; IDofs<NDofsMax; IDofs++) (*tExVec)[IDofs].val = ElmVecs(IElm, IDofs);

    //Perturb each of the DOF's
    for(int IDofs=0; IDofs<NDofsMax; IDofs++){
      (*tExVec)[IDofs].grad = 1.0;
      for(int IInteg=0; IInteg<NIntegs; IInteg++){
  //    ElmRess(IElm,IDofs) += (*tExVec)[IDofs].grad;
      }
      (*tExVec)[IDofs].grad = 0.0;
    }
  });

  //Get the residual vector and apply the essential BC's
  if(elem_restrict        != NULL) elem_restrict->MultTranspose(*EBlockResidual,y);
  if(ess_bcs_tdofs.Size() != 0)    y.SetSubVector(ess_bcs_tdofs,0.00);
};


/*****************************************\
!
!    Assemble and build the Jacobian
!               Operator
!
\*****************************************/
template<typename Number>
void tADNLForm<Number>::buildJacobian(const Vector & x) const
{
  //Get the element data vectors
  if(elem_restrict != NULL) elem_restrict->Mult(x,*EBlockVector);
  const auto ElmVecs = Reshape(EBlockVector->Read(), NDofsMax, nElms);

  mfem::forall_switch(use_dev, nElms, [=] MFEM_HOST_DEVICE (int IElm)
  {
    //Copy element vector in
    for(int IDofs=0; IDofs<NDofsMax; IDofs++) (*tEJVec)[IDofs].val.val = ElmVecs(IDofs, IElm);

    //Perturb each of the DOF's in Ith-row
    for(int IDofs=0; IDofs<NDofsMax; IDofs++){
      (*tEJVec)[IDofs].val.grad = 1.0;
      //Perturb each of the DOF's in Jth-column
      for(int JDofs=0; JDofs<NDofsMax; JDofs++){
        (*tEJVec)[JDofs].grad.val = 1.0;
        for(int IInteg=0; IInteg<NIntegs; IInteg++){
  //        elMats(IDofs, JDofs) += (*tEJVec).grad.grad;
        }
        (*tEJVec)[JDofs].grad.val = 0.0;
      }
      (*tEJVec)[IDofs].val.grad = 0.0;
    }
  });
};
