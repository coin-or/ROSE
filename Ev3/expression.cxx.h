/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        expression.h                                           *
* Source:      GNU C++                                                *
* Purpose:     symbolic expression (base classes and functionality)   *
* History:     010517 0.0 work started                                *
* License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
***********************************************************************/

#ifndef __EV3EXPRESSIONH__
#define __EV3EXPRESSIONH__

#define RCS5 "$Id: expression.h,v 1.25 2006/07/30 05:36:41 liberti Exp liberti $"

using namespace std;
#ifdef SUNWIN
#include <iostream>
#endif

#include <vector>
#include <string>
#include <map>

////////////// MEMCHECK debug ////////////////////
#ifdef MEMDEBUG
extern vector<pair<void*,pair<int,int*> > > memcheckdebug;
extern vector<pair<void*,pair<int,int*> > >::iterator memcheckdebugit;
extern pair<void*, pair<int, int*> > memcheckdebugpair;
extern int memcheckdebugcounter;
#endif
///////////// END MEMCHECK debug ////////////////

#include "common.h"
#include "tree.cxx"
#include "exceptions.h"
#include "fastexpression.h"

// algebraic expression operand
class Operand {

private:

protected:

  // one of the OperatorTypes above
  int oplabel;

  // 0 if no dependency, 1 if constant, 2 if coefficient, 3 if exponent
  int dependency;

  // if oplabel == CONST, the value of the constant
  double constant;

  // dependency for constants (added for MORON - see ../PROGNOTES)
  double* depconstant;

  // if oplabel == VAR, the index of the variable - should start from 1
  Int varindex;

  // if oplabel == VAR, the name of the variable
  string varname;

  // we allow multiplication for a constant coefficient in each Operand
  double coefficient;
  
  // dependency for coefficients (added for MORON - see ../PROGNOTES)
  double* depcoefficient;

  // we allow a real constant exponent in each Operand
  // THIS HAS MEANING ONLY IF operand IS A LEAF!!!
  double exponent;

  // dependency for exponents (added for MORON - see ../PROGNOTES)
  double* depexponent;

public:

  // constructors
  Operand();
  Operand(const double t);
  Operand(const Int t);
  Operand(const Int t, const bool isvar);
  // create a variable leaf and set coefficient
  Operand(const double c, const Int t, string vn);

  // copy constructor
  Operand(const Operand& t);

  // destructor
  ~Operand();

  // Operand class methods:

  // prints to a string
  string ToString(void) const;

  // get operator type
  int GetOpType(void) const;
  
  // get constant value - in CONSTs it multiplies by coefficient and
  // raises to exponent
  double GetValue(void) const;

  // just get the value, in all cases
  double GetSimpleValue(void) const;

  // get variable index
  Int GetVarIndex(void) const;

  // get variable name
  string GetVarName(void) const;

  // get the coefficient
  double GetCoeff(void) const;

  // get the exponent
  double GetExponent(void) const;

  // set operator type
  void SetOpType(const int t);
  
  // set constant value
  void SetValue(const double t);  

  // set variable index (start from 1 and add by steps of 1 when creating new
  // variables)
  void SetVarIndex(const Int t);

  // set variable name
  void SetVarName(const string vn);

  // set the exponent
  void SetExponent(const double expon);

  // set the coefficient
  void SetCoeff(const double coeff);

  // set constant dependencies (added for MORON - see ../PROGNOTES)
  void SetDependencyOnOperand(const int whichconstant, double** depvalue);

  // is operand a constant?
  bool IsConstant(void) const;

  // is operand a variable?
  bool IsVariable(void) const;

  // is operand a leaf node?
  bool IsLeaf(void) const;

  // is operand a zero constant?
  bool IsZero(void) const;
  
  // is operand a constant == v?
  bool HasValue(double v) const;

  // is operand a constant <= v?
  bool IsLessThan(double v) const;

  // is operand constant >= v?
  bool IsGreaterThan(double v) const;

  // set value = coefficient * value ^ exponent
  void ConsolidateValue(void);

  // enforce constant dependencies (added for MORON - see ../PROGNOTES)
  void EnforceDependencyOnOperand(void);

  // is operand this == operand t?
  bool operator == (const Operand& t);

  // substitute a variable with a constant
  void SubstituteVariableWithConstant(int varindex, double c);  

};

class BasicExpression;

typedef Pointer<BasicExpression> Expression;

class BasicExpression : public Operand, public Tree<BasicExpression> {

private:
  FastEvalTree* fetree;
  bool FastEvalAllocated;

public:

  // constructors

  // create empty
  BasicExpression();

  // create a constant leaf
  BasicExpression(const double t);

  // create a constant (integer-valued) leaf
  BasicExpression(const Int t);

  // create an (empty) operator or a variable leaf
  BasicExpression(const Int t, const bool isvar);
  
  // create a variable leaf and set coefficient
  BasicExpression(const double c, const Int t, string vn);

  // user-defined copy constructor with two options (to make a copy)
  BasicExpression(const Expression& t, const bool iscopy);

  // copy constructor
  BasicExpression(const BasicExpression& t);

  // destructor
  ~BasicExpression();

  // BasicExpression class methods:
  void Debug (void) const;

  // prints to a string
  string ToString(void) const;

  // output is a tree
  string PrintTree(int blanks, int tabs) const;

  // sets an expression to zero (deleting all existing subnodes)
  void Zero(void);

  // sets an expression to one (deleting all existing subnodes)
  void One(void);

  // is expression this == expression t?
  // (note that this half-replicates Tree::operator==,
  //  but I couldn't think of any other convenient way to fit in
  //  the operand data in == and still use the Tree's ==)
  bool IsEqualTo(const Expression& t) const;
  bool operator == (const BasicExpression& t) const;

  // other comparison operators
  bool IsEqualTo(const double t) const;

  // this returns true if args are equal up to top node coefficient
  bool IsEqualToNoCoeff(const Expression& t) const;

  // this returns true if args have equal operator structures and
  // maintain variable/constant assignment to leaf nodes
  // (a schema is an expression modulo leaves)
  bool IsEqualBySchema(const Expression& t) const;

  // this returns true if args have equal operator label
  bool IsEqualByOperator(const int theoplabel) const;

  // this returns the number of variables in the expression
  int NumberOfVariables(void) const;
  int NumberOfVariables(int& maxvi) const;

  // fast evaluation tree handling
  // creation (use nonrecursive version)
  void CreateFastEvalTree(void);
  void CreateFastEvalTreeRecursive(FastEvalTree* fet);
  // deletion (use nonrecursive version)
  void DeleteFastEvalTree(void);
  void DeleteFastEvalTreeRecursive(FastEvalTree* fet);
  
  // EVALUATION ROUTINES
  // Slow evaluators (for cases when just 1 evaluation is needed;
  // for repeated evaluations, use FastEval outside the class)
  // numeric evaluation: varvalues[i] contains the value of the
  // variable with varindex i, for all i < vsize.
  double Eval(double* varvalues, Int vsize) const;
  // in the following version, the varmap is used to map
  // var. indices in the expressions to indices of varvalues
  double Eval(double* varvalues, map<int,int>& varmap, Int vsize) const;
  // Fast evaluators (for repeated evaluations -- use nonrecursive version)
  double FastEval(double* varvalues, Int vsize);
  double FastEval(double* varvalues, map<int,int>& varmap, Int vsize);  
  double FastEval(FastEvalTree* fet, double* varvalues, Int vsize);
  double FastEval(FastEvalTree* fet, double* varvalues, 
		  map<int,int>& varmap, Int vsize);  
  double FastEvalRecursive(FastEvalTree* fet, double* varvalues, 
			   map<int,int>* varmap, Int vsize);

  // return the FastEvalTree pointer (build it if required)
  FastEvalTree* GetFastEvalTree(void);

  // whether expression depends on variable
  bool DependsOnVariable(Int vi) const;

  // whether expression depends on variable linearly
  // (0=depends nonlinearly, 1=depends linearly, 2=doesn't depend at all)
  int DependsLinearlyOnVariable(Int vi) const;

  // in a product, multiply all coeffs. of the operands, set result
  // as coeff of whole product, reset all operand coffs at 1; if
  // resulting global coeff is zero, delete all nodes and set 
  // this to zero constant.
  // don't do anything on other operators
  void ConsolidateProductCoeffs(void);

  // in a sum or product, if coeff of sum operand is not 1, distribute it
  // to the operands and set whole coeff to 1
  void DistributeCoeffOverSum(void);
  void DistributeCoeffOverProduct(void);

  // enforce constant dependencies (added for MORON - see ../PROGNOTES)
  // this only acts on the proper leaf nodes
  void EnforceDependency(void);

  // substitute a variable with a constant
  void VariableToConstant(int varindex, double c);

  // replace variable indexed v1 with variable indexed v2 (with varname vn)
  void ReplaceVariable(int v1, int v2, string vn);
  void ReplaceVariable(int v1, int v2, string vn, double c2);

  // replace with a variable the deepest node conforming to schema and
  // return replaced term or zero expression if no replacement occurs
  Expression ReplaceBySchema(int vi, string vn, Expression schema);
  // works on subnodes not on current node
  Expression ReplaceBySchemaRecursive(int vi, string vn,Expression schema);

  // replace with a variable the deepest node with given operator label
  // return replaced term or zero expression if no replacement occurs
  Expression ReplaceByOperator(int vi, string vn, int oplabel);  
  // works on subnodes not on current node
  Expression ReplaceByOperatorRecursive(int vi, string vn, int oplabel);  

  // replace all occurrences of subexpression needle with replace
  // return number of replacements
  int ReplaceSubexpression(Expression needle, Expression replace);

  // replace this with given expression (SIGSEGV risk, see implementation)
  void ReplaceWithExpression(Expression replace);

  // smith's standard form: takes a starting added variable index and
  // name and a vector of oplabels and one of schemata, performs replacements
  // on the this expression, returns sequence of defining constraints
  // for a defining constraint w_i = x o y, defcons[i] contains expression
  // x o y; caller must then use this to set constraint w - x o y = 0
  void SmithStandardForm(int defvarindex, string defvarname,
			 vector<int>& oplabels, vector<Expression>& schemata,
			 vector<Expression>& defcons);
  void SmithStandardFormRecursive(int defvarindex, string defvarname,
				  vector<int>& oplabels, 
				  vector<Expression>& schemata,
				  vector<Expression>& defcons);

  // prodbincont: replace products x*y where at least one var is binary
  // with w, and then adds exact linearization constraints for w
  // semantic difference of defcons with SmithStandardForm: here, 
  // defcons contains exactly 4 expressions lin_expr(w,x,y) 
  // to be set <= 0 by caller. Returns number of replacements.
  // addvarbounds contains map from added varindices to corresp var bounds
  int ProdBinCont(int defvarindex, string defvarname,
		  map<int,bool>& integrality,
		  map<int,double>& vlb, 
		  map<int,double>& vub,
		  map<int,pair<double,double> >& addvarbounds,
		  vector<Expression>& defcons);
  // returns number of replacements
  int ProdBinContRecursive(int defvarindex, string defvarname, 
			   map<int,bool>& integrality,
			   map<int,double>& vlb, 
			   map<int,double>& vub,
			   map<int,pair<double,double> >& addvarbounds,
			   vector<Expression>& defcons);

  // reset all names of variables having IDs between lid, uid to vn
  void ResetVarNames(string vn, int lid, int uid);
  
  // distribute products over sums - returns true if changed 
  // (re-call until false)
  bool DistributeProductsOverSums(void);

  // find variable indices in an expression
  void GetVarIndices(vector<int>& vidx);

  // return list of varIDs involved in a certain schema
  // e.g. f(x1*x2+x1*x3+x4*log(x5), x*y) = {1,2,3}
  void GetVarIndicesInSchema(vector<int>& vidx, Expression schema);

  // find the variable name corresponding to variable index vi
  string FindVariableName(int vi);

  // is this expression linear?
  bool IsLinear(void) const;

  // is this expression a quadratic product of variables? 
  // If yes, return the product type: PRODUCT, POWER or VAR
  bool IsQuadratic(int& prodtype) const;
  bool IsQuadratic(void) const;

  // return info about the linear part (assumes Simplify() has already been
  // called on this) - return false if expression has no linear part
  // by "linear part" we mean lin(x) in expr(x,y) = lin(x) + nonlin(y)
  // lincoeff[i] is the i-th nonzero coeff, linvi[i] is the corresponding
  // i-th varindex (starts from 1), c is the constant term to be added
  bool GetLinearInfo(vector<double>& lincoeff, 
		     vector<int>& linvi, vector<string>& linvn, double& c);

  // return info about the pure linear part (assumes Simplify() has
  // already been called on this) - much as above call but the "pure
  // linear part" is e.g. x+y in x+y+y^2
  bool GetPureLinearInfo(vector<double>& lincoeff, 
			 vector<int>& linvi, vector<string>& linvn, double& c);

  // get the linear part - x in x+y+y^2
  Expression GetLinearPart(void);

  // get the pure linar part - x+y in x+y+y^2
  Expression GetPureLinearPart(void);

  // get the nonlinear part - nonlin(y) in expr(x,y) = lin(x) + nonlin(y)
  Expression GetNonlinearPart(void);

  // get the purely nonlinear part - eg. y^2 in x+y+y^2
  Expression GetPureNonlinearPart(void);

  // return the additive constant of the expression 
  double GetConstantPart(void);

  // return the additive constant of the expression and remove it from the
  // expression itself
  double RemoveAdditiveConstant(void);

  // perform interval arithmetics on the expression;
  // Vlb[i] is the lower bound of varindex i, Vub is similar,
  // returned elb and eub hold values to expression interval
  void Interval(map<int,double> Vlb, map<int,double> Vub, 
		double& elb, double& eub);

  // this checks if expression is of the form
  // c * w_i  +/-  (unary_or_binary_operation) 
  // l and u are the constraint's lower and upper bounds
  // where the constraint is l < *this < u - it makes 
  // a difference for likening an inequality or equation
  // to Smith standard form
  bool IsSmithStandard(double l, double u);
  bool IsSmithStandard(void);

  // this checks if the set defined by the constraint 
  // l < *this < u is evidently convex 
  bool IsEvidentlyConvex(double l, double u);
  // this checks if expression *this is evidently convex 
  bool IsEvidentlyConvex(void);
  // this checks if expression *this is evidently concave 
  bool IsEvidentlyConcave(void);

  // this checks if expression is of the form
  // c * w_i  +/-  (operation for which we have a good convexifier) 
  bool IsOptStandard(void);

};

// All these functions contain tricks to simplify the operands. This
// means that both the operands may be changed, and indeed that the
// return value is often one of the changed operands.  To make sure
// you are not changing the operands, use the CopyOf function or the
// Copy() method in Pointer class; or use the equivalent operators below
// in order not to touch the argument expressions.

// add a link of b to a, return link to created expression
Expression SumLink(Expression a, Expression b);
// add a link of b to a with coeff -1, return link of a
Expression DifferenceLink(Expression a, Expression b);
// multiply a by a link of b, return link of a
Expression ProductLink(Expression a, Expression b);
// divide a by a link of b, return link of a
Expression FractionLink(Expression a, Expression b) throw(ErrDivideByZero);
// raise a to power b, return link of a
Expression PowerLink(Expression a, Expression b);
// change sign to coeff of a, return link of a
Expression MinusLink(Expression a);
// return link to log(a)
Expression LogLink(Expression a) throw(ErrNotPermitted);
// return link to exp(a)
Expression ExpLink(Expression a);
// other unary functions (of a):
Expression SinLink(Expression a);
Expression CosLink(Expression a);
Expression TanLink(Expression a);
Expression CotLink(Expression a) throw(ErrNotPermitted);
Expression SinhLink(Expression a);
Expression CoshLink(Expression a);
Expression TanhLink(Expression a);
Expression CothLink(Expression a) throw(ErrNotPermitted);
Expression SqrtLink(Expression a) throw(ErrNotPermitted);

// these are equivalent to the above but they don't change the arguments
Expression operator + (Expression a, Expression b);
Expression operator - (Expression a, Expression b);
Expression operator * (Expression a, Expression b);
Expression operator / (Expression a, Expression b) throw(ErrDivideByZero);
Expression operator ^ (Expression a, Expression b);
Expression operator - (Expression a);
Expression Log(Expression a) throw(ErrNotPermitted);
Expression Exp(Expression a);
Expression Sin(Expression a);
Expression Cos(Expression a);
Expression Tan(Expression a);
Expression Cot(Expression a) throw(ErrNotPermitted);
Expression Sinh(Expression a);
Expression Cosh(Expression a);
Expression Tanh(Expression a);
Expression Coth(Expression a) throw(ErrNotPermitted);
Expression Sqrt(Expression a) throw(ErrNotPermitted);

// symbolic differentiation of a w.r.t. variable with index vi, 
// return the created expression (a is not changed)
Expression Diff(const Expression& a, Int vi);
Expression DiffNoSimplify(const Expression& ac, Int vi);

// SIMPLIFICATIONS - all simplifications return true if they were
// effective or false if they weren't

// sin^2+cos^2 = 1 simplification 
bool TrigSimp(Expression a);

// generic simplification with modification of the expression
bool Simplify(Expression* a);
bool SimplifyRecursive(Expression* a);
bool SimplifyConstant(Expression* a);
bool DifferenceToSum(Expression* a);
bool CompactLinearPart(Expression* a);
bool CompactLinearPartRecursive(Expression* a);
bool CompactProducts(Expression* a);
bool ReorderNodes(Expression* a);

// destroys the whole tree and all nodes, be careful
void RecursiveDestroy(Expression* a);

// generic simplification on a copy of the expression
Expression SimplifyCopy(Expression* a, bool& ischanged);

  
#endif
