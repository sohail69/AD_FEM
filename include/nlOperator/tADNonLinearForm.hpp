#pragma once
#include "mfem.hpp"
#include "../UtilityObjects/macros.hpp"
#include "../UtilityObjects/lowLevelMFEM.hpp"
#include "../templatedMathObjs/dualNumber.hpp"
#include "../templatedMathObjs/tVector.hpp"
#include "../UtilityObjects/utilityFuncs.hpp"
#include "TCoeffInteg.hpp"
#include <vector>


//Linear algebra
//operators
#include "tInterpolator.hpp"
#include "tRestrictOperator.hpp"

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
  std::vector<TCoefficientIntegrator<Number>*>                                   EFuncCoeff;//Coeffs
  std::vector<std::function<dualSymNum<Number>(const tVector<dualSymNum<Number>> )>> dEfuncs;//dual-funcs
  std::vector<std::function<ddualSymNum<Number>(const tVector<ddualSymNum<Number>> )>> ddEfuncs;//ddual-funcs

  //Reference to block vector of element data
  Array<int> *btoffs_inp, *bvoffs_inp;
  Array<int> *btoffs_out, *bvoffs_out;
  mutable mfem::Vector  *EBlockVector, *EBlockResidual;
  mutable mfem::DenseMatrix elMats;

  //Used for directional derivatives Templated
  //dual number vector for Residual and Jacobian
  mutable tVector<dualSymNum<Number>>  *tERVec;
  mutable tVector<ddualSymNum<Number>> *tEJVec;
  mutable Operator *Jacobian_f=NULL;

  //Used for Isogeometric interpolator functions
  //at all integration points for a given integration
  //rule
//  Array<IntegrationRule> IntRules;

  //Problem sizing and restriction operator
  int NIntegsR=0, NIntegsJ=0;
  int nElms=0, NDofsMax=0, nElmDofs=0, NEQs=0;
  mfem::Operator *elem_restrict=NULL;
  tInterpolator IOp;


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
  void PrepareIntegrators();

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
  tERVec = new tVector<dualSymNum<Number>>(NDofsMax,mt);
  tEJVec = new tVector<ddualSymNum<Number>>(NDofsMax,mt);

  //////////////////////////
  ///Check sizes
  //////////////////////////
  std::cout <<  tERVec->size           << std::endl;
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
  delete tERVec, tEJVec, EBlockVector, EBlockResidual;
};


/*****************************************\
!
! Add energy functional coeffs
!
\*****************************************/
//Add a coefficient to the list
template<typename Number>
void tADNLForm<Number>::AddEnergyFuncCoeff(TCoefficientIntegrator<Number> *EFuncCoeff
                                         , std::vector<int> InpBlocks)
{


};

//Finalize the Integrator functions before assembly
template<typename Number>
void tADNLForm<Number>::PrepareIntegrators(){


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
  //with a restrict operator
  if(elem_restrict != NULL) elem_restrict->Mult(x,*EBlockVector);

  //Some stuff for devices
  const auto ElmVecs = Reshape(EBlockVector->Read(), NDofsMax, nElms);
  auto ElmRes = Reshape(EBlockResidual->ReadWrite(), NDofsMax, nElms);

  unsigned sum_NIpsNDofs=0;
  for(int IIntegs=0; IIntegs<; ) sum_NIpsNDofs
  tVector<dualSymNum<Number>> g_Xp, fg_Xp;

  //Partially assemble the sampled
  //variables used for calculating
  //the residuals
  mfem::forall_switch(use_dev, nElms*sum_NIpsNVars, [=] MFEM_HOST_DEVICE (int Ik)
  {
    //The recover the sub iterators
    unsigned IElm, IDof, Ip;
    InvIterator(Ik, IElm, IDof, Ip);

    //Interpolate the DOF's to
    //get the vars at the sample points 
    dualSymNum<Number> zero(0.00), dx(0.00,1.00);
    for(unsigned JDof=0; JDof<NDofMax; JDof++){
      g_Xp(Ik) = IOp.GetMat()(Ik,JDof)*(ElmVecs(IElm,JDof) + ((JDof==IDof)?zero:dx) + ((JDof==0)?zero:g_Xp(Ik));
    }

    fg_Xp(Ik) = func(g_Xp, Ik, IElm, IDof, Ip);
  });



  //Get the residual vector and apply the essential BC's
  if(elem_restrict        != NULL) elem_restrict->MultTranspose(*EBlockResidual,y);
  if(ess_bcs_tdofs.Size() != 0)    y.SetSubVector(ess_bcs_tdofs,0.00);
};

/*****************************************\
!
!     Assemble the Jacobian matrix
!
\*****************************************/
void buildJacobian(const Vector & x){};

};
