#pragma once
#include "mfem.hpp"
#include "../UtilityObjects/macros.hpp"
#include "../UtilityObjects/lowLevelMFEM.hpp"
#include "../templatedMathObjs/dualNumber.hpp"
#include "../templatedMathObjs/tVector.hpp"
#include "../UtilityObjects/utilityFuncs.hpp"
#include <vector>


//Linear algebra
//operators
#include "../templatedMathObjs/tMultiVarVector.hpp"
#include "tInterpolator.hpp"
#include "tRestrictOperator.hpp"

template<typename Num> using dualSymNum = dualNumber<Num,Num>;

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
! Basic note:
! TrueVars are Vars the have a discrete
! data analog, I.e. they have a contribution
! in the global vector and associated DOF's.
! Vars(ProblemVars) are a set of
! variables that can be derived by 
! interpolating the duscrete data-set of
! TrueVars
!
! example problem:
! (u + k.grad(u)) . v = f
! solved for u
! in this example u is a TrueVar, grad(u)
! is a Var/ProblemVar that can be derived
! by choice of interpolation of the discrete
! values of u
!
\*****************************************/
template<typename Number>
class tADNLForm : public Operator
{
private:
  //Reference to essential boundary conditions
  const std::vector<ParGridFunction*> & TrueVars;
  std::vector<mfem::Array<int>*>      ess_bcs_markers;
  mfem::Array<int>                    ess_bcs_tdofs;

  //Functions for evaluating the coefficients
  //at the at the integration points for residual
  //and the Jacobian
//TODO:Make a var map for Var blocks used in Jacbian forms
  std::vector<std::function<void(const mfem::Vector         & x
                               , const MFEMVarIterData<int> & Iter
                               , const mfem::Array<int>     & activeBlocks)>> Rfuncs; 


  std::vector<std::function<void(const Vector & x, const MFEMVarIterData<int> & Iter)>> Jfuncs; 

  //Reference to block vector of element data
  mutable mfem::Vector  *xE_Samp=NULL;
  mutable mfem::Vector  *EBlockVector=NULL, *EBlockResidual=NULL;
  mutable mfem::DenseMatrix elMats;

  //Used for directional derivatives Templated
  //dual number vector for Residual and Jacobian
  mutable mfem::Operator *Jacobian_f=NULL;

  //Restriction and Interpolation operators
  mutable mfem::Operator *elem_restrict=NULL;
  mutable tInterpolator *IOp=NULL;

  //Device and memory configs
  const bool             & use_dev;
  const mfem::Device     & device;
  const mfem::MemoryType & mt;

  //Problem sizing
  int nVars=0, sum_nIps_nVars=0;
  int nElms=0, nDofsMax=0, nElmDofs=0, nEQs=0;
  int OperatorSize(const std::vector<ParGridFunction*> & TrueVars_);

  //Iterators for MultiVarTensor data
  mutable bool VarIterUpdateFlag=false;
  VarIterData<int>     IO_VarIterator;
  MFEMVarIterData<int> MFEM_VarIterator;

public:
  //Constructor
  tADNLForm(const std::vector<ParGridFunction*> & TrueVars_, const mfem::Device & dev
          , const mfem::MemoryType & mt_, const bool & use_dev_);

  //Destructor
  ~tADNLForm();

  //Add a tensor variable to the list of sampled
  //variables which are sampled over every element
  //this variable has a parent True Var and an interpolator
  void AddTVar(const Var<int> & newVar, const mfem::Operator & InterpOp);

  //Prepare the operator before solving the
  //problem, does miscallaneous things such as:
  // ->Updates the MultiVarIterator
  // ->Updates the interpolator for the sampled Vars
  // ->Sizes the vectors needed for sampling the Vars
  void PrepareOperator() const;

  /// Assembles the residual form i.e.
  /// sums over all domain/bdr coefficients.
  virtual void Mult(const Vector & x, Vector & y) const;

  //Build the Jacobian Operator
  void buildJacobian(const Vector & x) const;

  //Returns a handle to the Jacobian
  mfem::Operator & GetGradient(const mfem::Vector &x) const override;
};

//An inverse iterator
//to 
template<typename uint>
void InvIterator(const uint & Ik, uint & IElm, uint & IDof, uint & Ip){

};

/*****************************************\
!
! Construct the non-linear form
!
\*****************************************/
//Size the operator
template<typename Number>
int tADNLForm<Number>::OperatorSize(const std::vector<ParGridFunction*> & TrueVars_)
{
  int Size=0;
  for(int I=0; I<TrueVars_.size(); I++) Size += TrueVars_[I]->ParFESpace()->TrueVSize();
  return Size;
};

//The constructor
template<typename Number>
tADNLForm<Number>::tADNLForm(const std::vector<ParGridFunction*> & TrueVars_, const mfem::Device & dev
                           , const mfem::MemoryType & mt_, const bool & use_dev_):
                             mfem::Operator(OperatorSize(TrueVars_),OperatorSize(TrueVars_))
                           , TrueVars(TrueVars_), device(dev), mt(mt_), use_dev(use_dev_)
{
  //////////////////////////
  ///Recover the object sizes
  //////////////////////////
  nElms = TrueVars[0]->ParFESpace()->GetMesh()->GetNE();
  nEQs  = OperatorSize(TrueVars_);
  for(int I=0; I<TrueVars.size(); I++){
    int NDofVar=0;
    for(int J=0; J<nElms; J++){
      DofTransformation doftrans;
      mfem::Array<int> dofs;
      TrueVars[I]->ParFESpace()->GetElementDofs(J,dofs,doftrans);
      NDofVar = std::max(dofs.Size(),NDofVar);
      nElmDofs += dofs.Size();
    }
    nDofsMax += NDofVar;
  }

  //////////////////////////
  ///Allocate the memory objects
  //////////////////////////
  elMats.SetSize(nDofsMax);
  EBlockVector   = new mfem::Vector(nElmDofs,mt_);
  EBlockResidual = new mfem::Vector(nElmDofs,mt_);

  //////////////////////////
  ///Check sizes
  //////////////////////////
  std::cout <<  EBlockVector->Size()   << std::endl;
  std::cout <<  EBlockResidual->Size() << std::endl;

  //////////////////////////
  ///Clear the Multi-Variate
  ///tensor sampled Vars
  // (for good measure)
  //////////////////////////
  clearIterator(IO_VarIterator);
};


/*****************************************\
!
! Destroy the Form
!
\*****************************************/
template<typename Number>
tADNLForm<Number>::~tADNLForm()
{
  clearIterator(IO_VarIterator);
  delete EBlockVector, EBlockResidual;
};

/*****************************************\
!
!  Adding in a variable to be sampled by
!  the assembly (FULL|ELEMENT|PARTIAL) loop
!  in the iterator class
!
\*****************************************/
template<typename Number>
void tADNLForm<Number>::AddTVar(const Var<int> & newVar, const mfem::Operator & InterpOp)
{
  //IOp->AddInterpolator(newVar.ParentTrueVar, InterpOp)
  AddVarIteratorDat(IO_VarIterator, newVar.TRank, newVar.sizes);
  VarIterUpdateFlag=true;
};

/*****************************************\
!
!  Preparing the operator for Mult
!  and GetGradient functions
!
\*****************************************/
template<typename Number>
void tADNLForm<Number>::PrepareOperator() const
{


}

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
  //Rebuild the MultiVarIterator,
  //Interpolator and sampler, if
  //the user has forgotten
  if(VarIterUpdateFlag) PrepareOperator();

  //Get the element data vectors
  //with a restrict operator
  if(elem_restrict != NULL) elem_restrict->Mult(x,*EBlockVector);

  //Some stuff for devices
  const auto ElmVecs = Reshape(EBlockVector->Read(), nDofsMax, nElms);
  auto ElmRes = Reshape(EBlockResidual->ReadWrite(), nDofsMax, nElms);
  mfem::Vector xtmp(nDofsMax,mt);

  //Partially assemble the sampled
  //variables used for calculating
  //the residuals
  mfem::forall_switch(use_dev, nElms*sum_nIps_nVars, [=] MFEM_HOST_DEVICE (int Ik)
  {
    //The recover the sub iterators
    unsigned IElm, IDof, Ip;
    InvIterator<unsigned>(Ik, IElm, IDof, Ip);


    //Interpolate the DOF's to get
    //the TrueVars at the sample points
    for(unsigned JDof=0; JDof<nDofsMax; JDof++){
 //     g_Xp(Ik) = IOp.GetMat()(Ik,JDof)*(ElmVecs(IElm,JDof) + ((JDof==0)?zero:g_Xp(Ik));
    }

  //for(unsigned ICoeff=0; ICoeff<Rfuncs.size(); ICoeff++){
  //    Rfuncs(xtmp, MFEM_VarIterator)
//  std::vector<std::function<void(const Vector & x, const MFEMVarIterData<int> & Iter)>> Rfuncs; 
// Rfuncs; 


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
template<typename Number>
void tADNLForm<Number>::buildJacobian(const Vector & x) const {};
