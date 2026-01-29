#pragma once
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
//#include "tCmath.hpp"


/*****************************************\
!
!  This is a parsing library for parsing
!  mathematical expressions found in the
!  tCmath library and basic operators
!
\*****************************************/
struct Token{
  std::string type;
  std::string value;
};

/*****************************************\
!
!  A simple lexer that is used for
!  tokenization and lexographic analysis
!  expressions
!
\*****************************************/
void lex(std::string & InputStr, std::vector<Token> & tokens){
  std::stringstream InpSS;
  InpSS << InputStr;
  std::string parseVal;
  while(std::getline(InpSS, parseVal, ' ')){
    Token nTk;
    nTk.value = parseVal;
    tokens.push_back(nTk);
  };
  std::cout << std::endl;
};

