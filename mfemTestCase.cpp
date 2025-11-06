#include <iostream>
#include <iomanip>
#include <cmath>
#include "mfem.hpp"
#include "../../MFEM_STUFF/mfem-4.7/build/include/mfem/general/forall.hpp"
#include "include/nlOperator/nlResidualForm.hpp"


// Mathematical objects
#include "include/templatedMathObjs/dualNumber.hpp"
#include "include/templatedMathObjs/complexNumber.hpp"
#include "include/templatedMathObjs/tVector.hpp"

// Mathematical Functions
#include "include/templatedMaths/tExpLogFuncs.hpp"
#include "include/templatedMaths/tPowFuncs.hpp"
#include "include/templatedMaths/tHypFuncs.hpp"
#include "include/templatedMaths/tErfGammaFuncs.hpp"

// Mathematical Functions

using ADdouble    = dualNumber<mfem::real_t,mfem::real_t>;
using JacADdouble = dualNumber<dualNumber<mfem::real_t,mfem::real_t>,dualNumber<mfem::real_t,mfem::real_t>>;

int main(){
  bool use_dev=true;
  const UINT64 nElms = 20, nDofs=8;
  const UINT64 n = nElms*nDofs;
  const char *device_config = "cpu";
  mfem::Device device(device_config);
  mfem::MemoryType mt = device.GetMemoryType();


  //resForm.Assemble(const Vector & x);
  mfem::Vector      EBlockVector(n,mt); //Standard Real_t vector that holds element block data
  tVector<ADdouble> tEVec(n,mt);        //Templated dual number vector for auto-diff
  residualForm<mfem::real_t> resForm(device, mt);
  resForm.Assemble(EBlockVector);


  auto d_Y = EBlockVector.Read(use_dev);
  auto d_X = tEVec.ReadWrite(use_dev);
  mfem::forall_switch(use_dev, n, [=] MFEM_HOST_DEVICE (int IElm)
  {
    for(UINT64 J=0; J<10; J++){ //Over all Integration points
      const ADdouble a(1.0,1.0), b(1.0,0.0);
      d_X[IElm] = exp<ADdouble>(a);
    }
  });

//  std::cout << tVec

  return 0;
};












