/*
** Name:       solver_Rsmith.h
** Author:     Leo Liberti 
** Source:     GNU C++
** Purpose:    Smith's standard form reformulation
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    080225 work started
*/

#ifndef _ROSESOLVERSMITHH
#define _ROSESOLVERSMITHH

#include "solver.h"

#define ROSEINFINITY 1e30

class SmithSolver : public virtual Solver {

 public:
  SmithSolver();
  ~SmithSolver();

  // data access
  void SetProblem(Problem* p); 
  std::string GetName(void) { return TheName; }
  Problem* GetProblem(void) { return InProb; }
  bool IsProblemInitialized(void) { return InitProb; }
  int GetOptimizationDirection(void) { return OptDir; }
  void SetOptimizationDirection(int theoptdir);
  bool CanSolve(int problemtype);

  // parameters
  Parameters GetParams(void) const { return ParameterBlob; }
  Parameters& GetParamsRef(void) { return ParameterBlob; }
  void ReplaceParams(Parameters& Blob) { ParameterBlob = Blob; }

  // initialize and solve
  void Initialize(bool force);
  int Solve(void) { return Solve(false); }
  int Solve(bool reinitialize);
  // get solution methods: if set to INFINITY, infeasible, if NaN, unsolved
  // get solution in maps ID -> value
  void GetSolution(std::map<int,double>& objfunval, std::map<int,double>& solution);
  // get solution in vectors ordered by local indices
  void GetSolutionLI(std::vector<double>& objfunval, std::vector<double>& solution);
  // as above but just 1 objective 
  void GetSolutionLI(double& objfunval, std::vector<double>& solution);
  
  // starting point
  void SetStartingPoint(int varID, double sp);
  void SetStartingPointLI(int localindex, double sp);
  double GetStartingPoint(int varID) const;
  double GetStartingPointLI(int localindex) const;
  // bounds
  void SetVariableLB(int varID, double LB);
  double GetVariableLB(int varID) const;
  void SetVariableUB(int varID, double UB);
  double GetVariableUB(int varID) const;
  void SetVariableLBLI(int localindex, double LB);
  double GetVariableLBLI(int localindex) const;
  void SetVariableUBLI(int localindex, double UB);
  double GetVariableUBLI(int localindex) const;

 private:
  // input optimization problem
  Problem* InProb;

  // name of solver
  std::string TheName;
  
  // is problem initialized?
  bool InitProb;

  // number of vars/constrs in the original problem
  int NumberOfVariables;
  int NumberOfConstraints;
  int NumberOfObjectives;

  // has this solver solved the problem yet?
  bool IsSolved;

  // is the most recent solution feasible in the problem?
  bool IsFeasible;

  // accuracy of solution
  double Epsilon;

  // name of local solver
  std::string LocalSolverName;
  
  // local solver
  Solver* LocalSolver;

  // maximum running time (disabled if = 0)
  double MaxRunningTime;

  // optimization direction (Minimization = 0, Maximization = 1)
  int OptDir;
  // optimization direction coefficient (1 for min, -1 for max)
  double OptDirCoeff;
  // true if opt dir has been set manually by the user (don't overwrite
  // with data from problem)
  bool ManualOptDir;

  // starting point
  std::vector<double> StartingPoint;

  // current solution
  std::vector<double> x;
  double f;
  std::vector<double> xstar;
  double fstar;

  // integrality
  std::vector<bool> integrality;

  // variable bounds
  std::vector<double> vlb;
  std::vector<double> vub;

  // optimal objective function value
  double OptimalObjVal;

  // parameters
  Parameters ParameterBlob;

  // be quiet
  bool Quiet;

  // name of main solver
  std::string MainSolver;
  
  // true if called from ampl
  bool AmplFlag;

};

#endif
