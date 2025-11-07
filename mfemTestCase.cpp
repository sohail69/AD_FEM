#include <iostream>
#include <iomanip>
#include <cmath>
#include "mfem.hpp"
#include "../../MFEM_STUFF/mfem-4.7/build/include/mfem/general/forall.hpp"
#include "include/nlOperator/nlForm.hpp"


// Mathematical objects
#include "include/templatedMathObjs/dualNumber.hpp"
#include "include/templatedMathObjs/complexNumber.hpp"
#include "include/templatedMathObjs/tVector.hpp"

// Mathematical Functions
//#include "include/templatedMaths/tCmath.hpp"


int main(){
  bool use_dev=false;
  const UINT64 nElms = 20, nDofs=8;
  const UINT64 n = nElms*nDofs;
  const char *device_config = "cpu";
  mfem::Device device(device_config);
  mfem::MemoryType mt = device.GetMemoryType();

  //Test case
  mfem::Vector x(n,mt), y(n,mt); 
  nlForm<mfem::real_t> nlProb(device, mt, use_dev);
  nlProb.buildJacobian(x);
  nlProb.Mult(x,y);

  return 0;
};












