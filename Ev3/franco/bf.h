/**********************************************************************
* Author:      Franco Raimondi                                        *
* Name:        bf.h                                                   *
* Source:      GNU C++                                                *
* Purpose:     Boolean expressions (base classes and functionality)   *
* History:     040609 work started                                    *
***********************************************************************/

#ifndef __BFH__
#define __BFH__

using namespace std;
#ifdef SUNWIN
#include <iostream>
#endif

#include <string>
#include <map>
#include "tree.cxx"
#include "exceptions.h"

#define NOVARIABLE -1
#define OPERANDSTRBUF 64
#define EXPRESSIONSTRBUF 2048

#define LARGE +1E10

// various operator types
enum OperatorType {  
  AND, OR, NOT,
  VAR, CONST,
  ERROR
};

// algebraic expression operand
class Operand {

private:

protected:

  // one of the OperatorTypes above
  int oplabel;

  // if oplabel == CONST, the value of the constant (true or false)
  bool constant;

  // if oplabel == VAR, the index of the variable - should start from 1
  int varindex;

  // if oplabel == VAR, the name of the variable
  string varname;

  // we allow multiplication for a constant coefficient in each Operand

public:

  // constructors
  Operand();
  Operand(const bool t); // True or False

  Operand(const int t, string vn); // A variable

  // copy constructor
  Operand(const Operand& t);

  // destructor
  ~Operand();

  // Operand class methods:

  // prints to a string
  string ToString(void) const;

  // get operator type
  int GetOpType(void) const;
  
  // get constant value
  bool GetValue(void) const;

  // get variable index
  int GetVarIndex(void) const;

  // get variable name
  string GetVarName(void) const;

  void SetOpType(const int t);
  
  // set constant value
  void SetValue(const bool t);  

  // set variable index (start from 1 and add by steps of 1 when creating new
  // variables)
  void SetVarIndex(const int t);

  // set variable name
  void SetVarName(const string vn);

  // is operand a constant?
  bool IsConstant(void) const;

  // is operand a variable?
  bool IsVariable(void) const;

  // is operand a leaf node?
  bool IsLeaf(void) const;

  // is operand a zero constant?
  bool IsZero(void) const;
  
  // is operand this == operand t?
  bool operator == (const Operand& t);

  // substitute a variable with a constant
  void SubstituteVariableWithConstant(long int varindex, bool c);  

};

class BasicBoolForm;

typedef Pointer<BasicBoolForm> BoolForm;

class BasicBoolForm : public Operand, public Tree<BasicBoolForm> {

private:

public:

  // constructors

  // create empty
  BasicBoolForm();

  // create a constant leaf
  BasicBoolForm(const bool t);

  // create a variable leaf
  BasicBoolForm(const int t, string vn);

// user-defined copy constructor with two options
  BasicBoolForm::BasicBoolForm(const BoolForm& t, const bool iscopy);

  // copy constructor
  BasicBoolForm(const BasicBoolForm& t);

  // destructor
  ~BasicBoolForm();

  // BasicBoolForm class methods:
  void Debug (void) const;

  // prints to a string
  string ToString(void) const;

  // output is a tree
  string PrintTree(int blanks, int tabs) const;

  // sets a formula to false (deleting all existing subnodes)
  void Zero(void);

  // sets a formula to true (deleting all existing subnodes)
  void One(void);

  // is expression this == expression t?
  // (note that this half-replicates Tree::operator==,
  //  but I couldn't think of any other convenient way to fit in
  //  the operand data in == and still use the Tree's ==)
  bool IsEqualTo(const BoolForm& t) const;
  bool operator == (const BasicBoolForm& t) const;

  // this returns the number of variables in the expression
  int NumberOfVariables(void) const;

  // substitute a variable with a constant
  void VariableToConstant(long int varindex, bool c);

  // replace variable indexed v1 with variable indexed v2 (with varname vn)
  void ReplaceVariable(long int v1, long int v2, string vn);

  // find the variable name corresponding to variable index vi
  string FindVariableName(long int vi);

};

BoolForm operator + (BoolForm a, BoolForm b);
BoolForm operator * (BoolForm a, BoolForm b);
BoolForm operator ! (BoolForm a);

BoolForm ToCNF( BoolForm a );

BoolForm PushNeg ( BoolForm a );

  
#endif
