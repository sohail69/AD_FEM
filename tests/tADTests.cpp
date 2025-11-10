#include <cmath>
#include <vector>
#include "../include/templatedMaths/tCmath.hpp"
#include "../include/templatedMathsObjs/dualNumber.hpp"


int main(){
  //
  // Testing my Trig functions
  // (outputs a dat file that
  // can be visualised with GNU-Plot)
  //
  std::ofstream OutFile("trigFuncTest.dat");
  const unsigned N=550;
  const double PI=3.14159265358979323846264338328;
  dualNumber<double,double> theta(0.0,1.0), zero(0.0,0.0);
  for(int I=0; I<N; I++){
    theta.val = double(5.0*I*PI/(N-1)) - PI;
    dualNumber<double,double> stheta = sin(theta);
    OutFile << theta.val << "   " << stheta.val          << "   " << stheta.grad 
                         << "   " << std::sin(theta.val) << "   " << std::cos(theta.val)  << std::endl;
  }
  OutFile.close();

  return 0;
};
