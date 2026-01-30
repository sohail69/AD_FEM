#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>


// Mathematical objects
/*#include "include/templatedMathObjs/dualNumber.hpp"
#include "include/templatedMathObjs/complexNumber.hpp"
#include "include/templatedMathObjs/tVector.hpp"

// Mathematical Functions
#include "include/templatedMaths/tCmath.hpp"*/
#include "include/templatedMaths/tParser/tCParser.hpp"


// Compile:
// mpic++ -m64 -std=c++17 -O2 -o main AD_ParseTest.cpp
int main(){
  std::string Iters    = "I J";
  std::string varSizes = "Dim1 Dim2";
  std::string Vars     = "a U[Dim2] V[Dim1] D[Dim1,Dim2]";
  std::string expr     = " a + D[J,I]*U[I]*V[J]";

  BTreeNode parseTree;
  std::vector<double> inpData;
  std::vector<Token> tokenizedData;
  std::map<unsigned,std::function<double(double)>> funcMap;
  auto lmbdaFunc = tensorParse<double>(Iters, varSizes, Vars, expr);
  return 0;
};

/*
//  std::cout << lmbdaFunc(inpData) << std::endl;
    std::cout << std::setw(15) << VarTokens[I].type
              << std::setw(15) << VarTokens[I].size
              << std::setw(15) << VarTokens[I].value << std::endl;
*/
