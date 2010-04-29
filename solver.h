/*
** Name:       solver.h
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    "Common" problem solver header file
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    050217 work started
*/

#ifndef _ROSESOLVERH
#define _ROSESOLVERH

#include "problem.h"
#include <vector>

#define NOTIMPL std::cerr << "method not implemented\n"; exit(11)

enum ProblemType { ROSE_LP, ROSE_MILP, ROSE_NLP, ROSE_MINLP };

class UserCut {
 public:
  UserCut() {}
  UserCut(Expression e, double theLB, double theUB);
  Expression Function;
  FastEvalTree* FunctionFET;
  Expression NonlinearPart;
  FastEvalTree* NonlinearPartFET;
  double LB;
  double UB;
  // is cut linear?
  bool IsLinear;
  // we also record first and second symbolic derivatives (if needed)
  std::vector<Expression> Diff;
  std::vector<FastEvalTree*> DiffFET;
  std::vector<std::vector<Expression> > Diff2;
  std::vector<std::vector<FastEvalTree* > > Diff2FET;
};

class UserLinearCut {
 public:
  UserLinearCut();
  UserLinearCut(const UserLinearCut& uc);
  UserLinearCut(std::vector<std::pair<int, double> >& lininfo, 
		double theLB, double theUB);
  UserLinearCut(int* varindices, double* coeffs, int asize, 
		double theLB, double theUB);
  ~UserLinearCut();
  UserLinearCut& operator=(const UserLinearCut& uc);
  double LB;
  double UB;
  int Nonzeroes;
  int* VarIndices;
  double* Coeffs;
  int* memcount;
};

class Solver {

 public:
  Solver() { }
  ~Solver() { }

  // general methods
  virtual int Solve(void) { NOTIMPL; }
  virtual void Initialize(bool force = false) { NOTIMPL; }
  virtual void SetProblem(Problem* p) { NOTIMPL; }
  virtual Problem* GetProblem(void) { NOTIMPL; return NULL; }
  virtual bool IsProblemInitialized(void) { NOTIMPL; return false; }
  virtual int GetOptimizationDirection(void) { NOTIMPL; return 0; }
  virtual void SetOptimizationDirection(int theoptdir) { NOTIMPL; }
  virtual void GetSolution(std::map<int,double>& objfunval, 
			   std::map<int,double>& solution) { NOTIMPL; }
  virtual void GetSolutionLI(std::vector<double>& objfunval, 
			     std::vector<double>& solution) { NOTIMPL; }
  virtual void GetSolutionLI(double& objfunval, std::vector<double>& solution) {
    NOTIMPL;
  }
  virtual std::string GetName(void) { NOTIMPL; }
  virtual bool CanSolve(int problemtype) { NOTIMPL; }

  // parameters
  virtual Parameters GetParams(void) const { NOTIMPL; }
  virtual Parameters& GetParamsRef(void) { NOTIMPL; }
  virtual void ReplaceParams(Parameters& Blob) { NOTIMPL; }

  // solver-specific cuts
  virtual void SetMaxNumberOfCuts(int mncuts) { NOTIMPL; }
  virtual int GetMaxNumberOfCuts(void) const { NOTIMPL; return 0; }
  virtual int GetNumberOfCuts(void) const { NOTIMPL; return 0; }
  virtual int AddCut(Expression e, double LB, double UB) { NOTIMPL; return 0; }
  virtual int AddCut(Expression e, double LB, double UB, 
		     bool diff2) { NOTIMPL; return 0; }
  virtual int AddCut(std::vector<std::pair<int, double> >& cutinfo,
		     double LB, double UB) { NOTIMPL; return 0; }
  virtual double EvalCut(int cutID, double* xval) { NOTIMPL; return 0; }
  virtual double EvalNLCut(int cutID, double* xval) { NOTIMPL; return 0; }
  virtual double EvalCutDiff(int cutID, int varID, double* xval) {
    NOTIMPL; return 0;
  }
  virtual double EvalCutDiffNoConstant(int cutID, int varID, double* xval) {
    NOTIMPL; return 0;
  }
  virtual double EvalCutDiff2(int cutID, int varID1, int varID2, double* xval){
    NOTIMPL; return 0;
  }
  virtual bool IsCutLinear(int cutID) { NOTIMPL; return false; }
  virtual void EnableCut(int cutID) { NOTIMPL; }
  virtual void DisableCut(int cutID) { NOTIMPL; }
  virtual void SetCutLB(int cutID, double LB) { NOTIMPL; }
  virtual double GetCutLB(int cutID) { NOTIMPL; return 0; }
  virtual void SetCutUB(int cutID, double UB) { NOTIMPL; }
  virtual double GetCutUB(int cutID) { NOTIMPL; return 0; }

  // modify bounds
  virtual void SetVariableLB(int varID, double LB) { NOTIMPL; }
  virtual double GetVariableLB(int varID) const { NOTIMPL; return 0;}
  virtual void SetVariableUB(int varID, double UB) { NOTIMPL; }
  virtual double GetVariableUB(int varID) const { NOTIMPL; return 0;}
  virtual void SetConstraintLB(int constrID, double LB) { NOTIMPL; }
  virtual double GetConstraintLB(int constrID) const { NOTIMPL; return  0;}
  virtual void SetConstraintUB(int constrID, double UB) { NOTIMPL; }
  virtual double GetConstraintUB(int constrID) const { NOTIMPL; return 0;}
  virtual void SetVariableLBLI(int localindex, double LB) { NOTIMPL; }
  virtual double GetVariableLBLI(int localindex) const { NOTIMPL; return 0;}
  virtual void SetVariableUBLI(int localindex, double UB) { NOTIMPL; }
  virtual double GetVariableUBLI(int localindex) const { NOTIMPL; return 0;}
  virtual void SetConstraintLBLI(int localindex, double LB) { NOTIMPL; }
  virtual double GetConstraintLBLI(int localindex) const { NOTIMPL; return 0;}
  virtual void SetConstraintUBLI(int localindex, double UB) { NOTIMPL; }
  virtual double GetConstraintUBLI(int localindex) const { NOTIMPL; return 0;}
  
  // starting point
  virtual void SetStartingPoint(int varID, double sp) { NOTIMPL; }
  virtual void SetStartingPointLI(int localindex, double sp) { NOTIMPL; }
  virtual double GetStartingPoint(int varID) const { NOTIMPL; return 0;}
  virtual double GetStartingPointLI(int localindex) const {
    NOTIMPL; return 0; 
  }

  // Lagrange multipliers of constraints and added cuts
  virtual double GetConstraintLagrangeMultiplier(int constrID) const 
    { NOTIMPL; return 0; }
  virtual double GetConstraintLagrangeMultiplierLI(int localindex) const
    { NOTIMPL; return 0; }
  virtual double GetCutLagrangeMultiplier(int cutID) const
    { NOTIMPL; return 0; }
  virtual double GetBoundLagrangeMultiplier(int varID) const
    { NOTIMPL; return 0; }
  virtual double GetBoundLagrangeMultiplierLI(int localindex) const
    { NOTIMPL; return 0; }

  // variable information
  virtual bool IsBasic(int varID) { NOTIMPL; return false; }
  virtual bool IsBasicLI(int localindex) { NOTIMPL; return false; }

  // simplex tableau
  virtual double GetSimplexTableauRow(int rowindex,
				      std::vector<std::pair<int,double> >& t) 
    { NOTIMPL; return 0; }

  // reformulator solver methods
  virtual int GetNumberOfReformulations(void) const { NOTIMPL; return 0; }
  virtual Problem* GetReformulation(int ridx) { NOTIMPL; return NULL; }
};


#endif
