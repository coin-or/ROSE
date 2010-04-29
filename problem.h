/*
** Name:       problem.h
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    header file defining the Problem class
** License:    (C) Leo Liberti, all rights reserved. Code published under the
               Common Public License.
** Note:       Problem entities are: variables, objectives, constraints;
               these are indexed in two ways: by their unique ID and by
               a "local index" which is relative to the number of such
               entities; e.g. a problem might have four variables with
               varIDs 3,5,10,11 whose local indices are 1,2,3,4 (because
               there are four of them); if varID 5 is deleted, the problem
               then has 3 variables whose IDs are 3,10,11 and whose local
               indices are 1,2,3
** History:    050215 work started
*/

#ifndef _ROSEPROBLEMH
#define _ROSEPROBLEMH

#define EPSILONTOLERANCE 1e-30
#define ROSEINFINITY 1e30
#define FIRSTOBJ 1

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <unistd.h>
#include <sys/time.h>
#include "Ev3/expression.h"
#include "Ev3/parser.h"
#include "parameters.h"
#include "utils.h"

// for now, keep this in, but in time we should get rid of this
using namespace std;

// ampl includes
// to add an AMPL option, also change problem.cxx, rose.cxx
namespace Ampl {
#include "asl.h"
#include "nlp.h"
#include "getstub.h"
#include "r_opn.hd" // for N_OPS
#include "opcode.hd"
  extern Option_Info* TheOInfo;
  extern int* TheConvexifierOutAmpl;
  extern int* TheGLPKSolutionMethod;
  extern double* TheLimitedBranchArbitrary;
  extern int* TheLimitedBranchLocalSolver;
  extern int* TheLimitedBranchMaxTime;
  extern int* TheLocalBranchKmax;
  extern int* TheLocalBranchWarmUp;
  extern int* TheMainSolver;
  extern int* TheNoSimplifier;
  extern int* TheQuiet;
  extern int* TheRelaxInteger;
  extern int* TheRQuarticConvexType;
  extern int* TheSnopt6MajorIterations;
  extern int* TheSnopt6MinorIterations;
  extern int* TheSymmgroupOutType;
  extern double* TheVNSEpsilon;
  extern int* TheVNSKmax;
  extern int* TheVNSKmin;
  extern int* TheVNSLocalSolver;
  extern int* TheVNSMaxTime;
  extern int* TheVNSNeighbourhood;
  extern int* TheVNSSamples;
  extern int* TheVNSWarmUp;
  extern int* TheDisaggrK;
  extern ASL* asl;
};

// ampl preliminary stuff
#undef exit
#define CHR (char*)  // for suppressing "String literal to char*" warnings
#define R_OPS ((ASL_fg*)asl)->I.r_ops_
#define OBJ_DE ((ASL_fg*)asl)->I.obj_de_
#define CON_DE ((ASL_fg*)asl)->I.con_de_

// for the optimization direction of the problem
enum {Minimization, Maximization};

// for the Ev3 parser
enum { NoSectionType, VariableSectionType, ObjfunSectionType,
       ConstraintSectionType, StartingpointSectionType,
       ParameterSectionType };

struct Variable {
  // variables consist of an ID, a name, integrality and a range
  int ID;
  std::string Name;
  double LB;
  double UB;
  bool IsIntegral;
  // flag to signal that variable cannot be deleted from problem
  bool Persistent;
  // whether variable is a constant, ie. not a decision variable
  bool IsConstant;
  // we also record optimal variable values
  double Optimum;
};

struct Objective {
  // objective functions consist of an ID, a name and an expression
  int ID;
  std::string Name;
  Expression Function;
  FastEvalTree* FunctionFET;
  Expression NonlinearPart;
  FastEvalTree* NonlinearPartFET;
  // optimization direction (Minimization=0, Maximization=1)
  int OptDir;
  // we also record first and second symbolic derivatives (if needed)
  std::vector<Expression> Diff;
  std::vector<FastEvalTree*> DiffFET;
  std::vector<std::vector<Expression> > Diff2;
  std::vector<std::vector<FastEvalTree*> > Diff2FET;
  // we also record the optimal objective function values
  double Optimum;
};

struct Constraint {
  // constraints consist of an ID, a name, an expression, and bounds
  int ID;
  std::string Name;
  Expression Function;
  FastEvalTree* FunctionFET;
  Expression NonlinearPart;
  FastEvalTree* NonlinearPartFET;
  double LB;
  double UB;
  // we also record first and second symbolic derivatives (if needed)
  std::vector<Expression> Diff;
  std::vector<FastEvalTree*> DiffFET;
  std::vector<std::vector<Expression> > Diff2;
  std::vector<std::vector<FastEvalTree*> > Diff2FET;
};

// predicates for finding variables, objectives, constraints given an ID
template<class T> class IsIDEqual {
 public:
  IsIDEqual() { }
  IsIDEqual(int theID) : myID(theID) { }
  ~IsIDEqual() { }
  void SetID(int theID) { myID = theID; }
  int GetID(void) const { return myID; }
  bool operator() (T& theEntity) {
    return myID == theEntity->ID;
  }
 private:
  int myID;
};

class Problem {

 public:

  // constructor, destructor
  Problem();
  Problem(bool thenosimplifierflag);
  ~Problem();

  // get general problem info
  std::string GetName(void) { return Name; }
  void SetName(std::string thename) { Name = thename; }
  bool IsProblemContinuous(void) { return IsProbContinuous; }
  bool IsProblemLinear(void) { return IsProbLinear; }
  Problem* GetParent(void) { return Parent; }
  int GetNumberOfChildren(void) { return Children.size(); }
  Problem* GetChild(int childindex);
  std::string GetFormulationName(void) { return FormulationName; }
  int GetOptimizationDirection(int objID);
  void SetOptimizationDirection(int objID, int theOptDir);
  int GetOptimizationDirectionLI(int localindex);
  void SetOptimizationDirectionLI(int localindex, int theOptDir);
  bool HasDeleted(void) { return DeletionFlag; }
  bool IsSolved(void) { return Solved; }
  void SetSolved(bool issolved) { Solved = issolved; }
  int IsFeasible(void) { return Feasible; }
  void SetFeasible(int isfeasible) { Feasible = isfeasible; }

  // parameters
  Parameters GetParams(void) const { return ParameterBlob; }
  Parameters& GetParamsRef(void) { return ParameterBlob; }
  void ReplaceParams(Parameters& Blob) { ParameterBlob = Blob; }

  // get problem sizes
  int GetNumberOfVariables(void) { return NumberOfVariables; }
  int GetNumberOfObjectives(void) { return NumberOfObjectives;}
  int GetNumberOfConstraints(void) { return NumberOfConstraints; }
  int GetNumberOfIntegerVariables(void);

  // get problem entities
  Variable* GetVariable(int varID);
  Variable* GetVariableLI(int localindex);
  Objective* GetObjective(int objID);
  Objective* GetObjectiveLI(int localindex);
  Constraint* GetConstraint(int constrID);
  Constraint* GetConstraintLI(int localindex);
  int GetVariableID(int localindex);  // inverse maps in Get---LocalIndex()
  int GetVarLocalIndex(int varID);
  int GetObjectiveID(int localindex);
  int GetObjLocalIndex(int objID);
  int GetConstraintID(int localindex);
  int GetConstrLocalIndex(int constrID);

  // get/set numerical info (with varID / with localindex)
  // optimal variable values
  double GetOptimalVariableValue(int varID);
  double GetOptimalVariableValueLI(int localindex);
  void SetOptimalVariableValue(int varID, double value);
  void SetOptimalVariableValueLI(int localindex, double value);
  // variable bounds
  void SetVariableLBValue(int varID, double value);
  void SetVariableLBValueLI(int localindex, double value);
  void SetVariableUBValue(int varID, double value);
  void SetVariableUBValueLI(int localindex, double value);
  //constant variables
  void SetConstantVariable(int varID, bool constant);
  void SetConstantVariableLI(int localindex, bool constant);

  // current variable values
  double GetCurrentVariableValue(int varID);
  double GetCurrentVariableValueLI(int localindex);
  void SetCurrentVariableValue(int varID, double value);
  void SetCurrentVariableValueLI(int localindex, double value);
  // test feasibility of a single constraint (indexed cID)
  bool TestConstraintsFeasibility(int cID, double tolerance, double& disc);
  double TestConstraintsInfeasibilityQuantity(int cID, double tolerance, double& disc);
  // test feasibility of current variables in problem constraints
  bool TestConstraintsFeasibility(double tolerance, double& discrepancy);
  // test feasibility of a single variable (indexed vID)
  bool TestVariablesFeasibility(int vID, double tolerance, double& disc);
  // test feasibility of current variables in problem ranges & integrality cnst
  bool TestVariablesFeasibility(double tolerance, double& discrepancy);
  // starting point
  double GetStartingPoint(int varID);
  double GetStartingPointLI(int localindex);
  // optimal obj fun values
  double GetOptimalObjectiveValue(int objID);
  void SetOptimalObjectiveValue(int objID, double value);
  // get solution methods: if set to INFINITY, infeasible, if NaN, unsolved
  // get solution in maps ID -> value
  void GetSolution(std::map<int,double>& objfunval, std::map<int,double>& solution);
  // get solution in vectors ordered by local indices
  void GetSolutionLI(std::vector<double>& objfunval, std::vector<double>& solution);

  // get additive constant of 1st objective
  double GetObj1AdditiveConstant(void) { return Obj1AdditiveConstant; }

  // evaluation subroutines
  double EvalObj(int objID);
  double EvalNLObj(int objID);
  double EvalObjDiff(int objID, int varID);
  double EvalObjDiffNoConstant(int objID, int varID);
  double EvalObjDiff2(int objID, int varID1, int varID2);
  double EvalConstr(int constrID);
  double EvalNLConstr(int constrID);
  double EvalConstrDiff(int constrID, int varID);
  double EvalConstrDiffNoConstant(int constrID, int varID);
  double EvalConstrDiff2(int constrID, int varID1, int varID2);
  bool IsObjConstant(int objID);
  bool IsObjDiffConstant(int objID, int varID);
  bool IsObjDiff2Constant(int objID, int varID1, int varID2);
  bool IsConstrConstant(int constrID);
  bool IsConstrDiffConstant(int constrID, int varID);
  bool IsConstrDiff2(int constrID, int varID1, int varID2);
  bool IsConstrActive(int constrID, double tolerance, int& LowerUpper);

  // read in problem from .ev3 file
  void Parse(char* thefilename);
  void ParseEv3String(char* element, int status, int lineno);

  // read in problem from AMPL
  Ampl::ASL* ParseAMPL(char** argv, int argc);
  Expression ParseAMPLExpressionTree(Ampl::expr* e);
  bool IsAMPLExpressionZero(Ampl::expr* e);

  // symbolic simplifier
  void Simplifier(bool theamplflag);
  void DebugSymbolicDiffs(void); // for debugging purposes only

  // return next available entity ID
  // (does not increase counters nor add the entity)
  int NewVariableID(void);
  int NewObjectiveID(void);
  int NewConstraintID(void);

  // entity creation (calls NewEntity() automatically within)
  void AddVariable(std::string& myName, bool myIntegral, bool myPersistent,
		   double LB, double UB, double myOptimum);

  void AddVariable(std::string& myName, bool myIntegral, bool myPersistent,
		   bool myConstant, double myLB, double myUB, double myOptimum);

  void AddObjective(std::string& myName, Expression myExpr, int myOptDir,
		    double myOptimum);
  void AddConstraint(std::string& myName, Expression myExpr, double LB,
		     double UB);

  // entity deletion
  // caution: deleting a variable from a problem does not delete it from
  // all the expressions occurring therein (i.e. this is not a projection)
  void DeleteVariable(int varID);
  void DeleteObjective(int objID);
  void DeleteConstraint(int constrID);

 protected:
  // problem name
  std::string Name;
  // flag signalling if problem only has continuous variables
  bool IsProbContinuous;
  // flag signalling if problem only has linear forms in objs and constrs
  bool IsProbLinear;
  // link to parent problem
  Problem* Parent;
  // link to children problem (reformulations)
  std::vector<Problem*> Children;
  // name of the form this problem is in
  std::string FormulationName;
  // true if the problem has already been solved by some numerical solver
  bool Solved;
  // 1 if feasible, -1 if infeasible, 0 if unknown
  int Feasible;

  // total number of variables
  int NumberOfVariables;
  // total number of objective functions
  int NumberOfObjectives;
  // total number of constraints
  int NumberOfConstraints;
  // variables
  std::vector<Variable*> Var;
  // objectives
  std::vector<Objective*> Obj;
  // constraints
  std::vector<Constraint*> Constr;
  // contains the next free variable ID
  int NextFreeVarID;
  // contains the next free objective ID
  int NextFreeObjFunID;
  // contains the next free constraint ID
  int NextFreeConstrID;

  // variable values in the form of an array of doubles
  double* x;
  // starting point
  std::vector<double> StartingPoint;
  // additive constant of objective function 1
  double Obj1AdditiveConstant;

  // maps between global and local indexing
  std::map<int,int> VarIDToLocalIndex;
  std::map<int,int> ObjIDToLocalIndex;
  std::map<int,int> ConstrIDToLocalIndex;

  // don't simplify
  bool NoSimplifier;

  // parameters
  Parameters ParameterBlob;

 private:
  // set if reformulation process has deleted variables or constraints
  bool DeletionFlag;

  // the expression parser used in Parse()
  ExpressionParser exprparser;

  // have diffs / diffs2 already been computed? if no, compute them
  bool ComputedSymbolicDiff;
  void ComputeSymbolicDiffs(void);
  bool ComputedSymbolicDiff2;
  void ComputeSymbolicDiffs2(void);
};

std::ostream& operator<<(std::ostream& s, Problem& gp);

#endif
