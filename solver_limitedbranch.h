/*
** Name:       solver_limitedbranch.h
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    solver limitedbranch header 
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    051110 work started
*/

#ifndef _ROSESOLVERLIMITEDBRANCHH
#define _ROSESOLVERLIMITEDBRANCHH

#include "solver.h"

#define ROSEINFINITY 1e30

class LimitedBranchSolver : public virtual Solver {

 public:
  LimitedBranchSolver();
  ~LimitedBranchSolver();

  // data access
  void SetProblem(Problem* p); 
  string GetName(void) { return TheName; }
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
  void GetSolution(map<int,double>& objfunval, map<int,double>& solution);
  // get solution in vectors ordered by local indices
  void GetSolutionLI(vector<double>& objfunval, vector<double>& solution);
  // as above but just 1 objective 
  void GetSolutionLI(double& objfunval, vector<double>& solution);
  
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

  // user cuts: these methods only change the problem data within 
  // the solver, not within the problem object itself
  void SetMaxNumberOfCuts(int mncuts); // only possible if uninitialized
  int GetMaxNumberOfCuts(void) const { return MaxNumberOfCuts; }
  int GetNumberOfCuts(void) const { return NumberOfCuts; }
  // adds a user linear cut, returns the cut ID; cuts are enabled by default
  int AddCut(std::vector<std::pair<int, double> >& cutinfo, 
	     double LB, double UB); 
  // adds a user cut, returns the cut ID; cuts are enabled by default
  int AddCut(Expression e, double LB, double UB); 
  // handle cuts
  void EnableCut(int cutID); // sets LB and UB to original values
  void DisableCut(int cutID); // sets LB to -inf and UB to inf
  void SetCutLB(int cutID, double LB);
  double GetCutLB(int cutID);
  void SetCutUB(int cutID, double UB);
  double GetCutUB(int cutID);

 protected:
  int Branch(int localvarlindex, 
	     std::vector<double> vlb, std::vector<double> vub);

 private:
  // input optimization problem
  Problem* InProb;

  // name of solver
  string TheName;
  
  // is problem initialized?
  bool InitProb;

  // number of vars/constrs in the original problem
  int NumberOfVariables;
  int NumberOfConstraints;

  // has this solver solved the problem yet?
  bool IsSolved;

  // is the most recent solution feasible in the problem?
  bool IsFeasible;

  // accuracy of solution
  double Epsilon;

  // name of local solver
  string LocalSolverName;
  
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
  vector<double> StartingPoint;

  // current solution
  vector<double> x;
  double f;
  vector<double> xstar;
  double fstar;

  // integrality
  vector<bool> integrality;

  // variable bounds
  vector<double> vlb;
  vector<double> vub;

  // optimal objective function value
  double OptimalObjVal;

  // parameters
  Parameters ParameterBlob;

  // be quiet
  bool Quiet;

  // relax integer variables?
  bool Relaxed;

  // timer
  double tstart[2], tend[2];
  double starttime, stoptime, usertimelap;

  // probability of making an arbitrary choice about left/right
  // branching precedence even though the objective function
  // usefully depends on the variable being branched on
  double ArbitraryChoice;

  // number of added cuts
  int NumberOfCuts;
  
  // maximum number of allowed added cuts (defaults to 20)
  int MaxNumberOfCuts;

  // do we need cuts?
  bool NeedsCuts;

  // current number of integer feasible solutions during branching
  int MINLPSols;

  // name of main solver
  string MainSolver;
  
  // whether rose called from ampl or not
  bool AmplFlag;

};

#endif
