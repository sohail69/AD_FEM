#include <fstream>
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
#include "include/templatedMaths/tCmath.hpp"


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


  //
  // Testing my Trig functions
  //
  std::ofstream OutFile("trigFuncTest.dat");
  const unsigned N=550;
  const double PI=3.14159265358979323846264338328;
  dualSymNum<mfem::real_t>  theta(0.0,1.0), zero(0.0,0.0);
  for(int I=0; I<N; I++){
    theta.val = double(25.0*I*PI/(N-1)) - PI;
    dualSymNum<mfem::real_t> stheta = sin(theta);
    OutFile << theta.val << "   " << stheta.val          << "   " << stheta.grad 
                         << "   " << std::sin(theta.val) << "   " << std::cos(theta.val)  << std::endl;
  }
  OutFile.close();

  return 0;
};
