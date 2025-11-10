#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "../include/templatedMaths/tCmath.hpp"

int main(){
  //
  // Testing my Trig functions
  // (outputs a dat file that
  // can be visualised with GNU-Plot)
  //
  std::ofstream OutFile("trigFuncTest.dat");
  const unsigned N=550;
  const double PI=3.14159265358979323846264338328;
  double theta(0.0);
  for(int I=0; I<N; I++){
    theta = double(5.0*I*PI/(N-1)) - PI;
    OutFile << theta.val << "   " << sin(theta)      << "   " << cos(theta) 
                         << "   " << std::sin(theta) << "   " << std::cos(theta)  << std::endl;
  }
  OutFile.close();


  return 0;
};
