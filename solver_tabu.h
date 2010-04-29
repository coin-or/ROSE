/*
** Name:       solver_tabu.h
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    GO tabu search solver header file (LETSGO)
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    050429 work started
*/

#ifndef _ROSESOLVERTABUH
#define _ROSESOLVERTABUH

#include "solver.h"

#define ROSEINFINITY 1e30

class TabuSolver : public virtual Solver {

 public:
  TabuSolver();
  ~TabuSolver();

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
  // as above but with just 1 objective
  void GetSolutionLI(double& objfunval, vector<double>& solution);
  
  // starting point
  void SetStartingPoint(int varID, double sp);
  void SetStartingPointLI(int localindex, double sp);
  double GetStartingPoint(int varID) const;
  double GetStartingPointLI(int localindex) const;

  // places the same random vector in x[] and localsolver's starting point
  void RandomLocalSolverStartingPoint(void); 

  // disturbs the starting point
  void DisturbLocalSolverStartingPoint(void);

  // finds a random unit direction
  void FindRandomUnitDirection(std::vector<double>& d);

 private:
  // solver name
  string TheName;

  // input optimization problem
  Problem* InProb;
  
  // is problem initialized?
  bool InitProb;

  // number of vars/constrs in the original problem
  int NumberOfVariables;
  int NumberOfConstraints;

  // has this solver solved the problem yet?
  bool IsSolved;

  // epsilon tolerance
  double Epsilon;

  // optimization direction (Minimization = 0, Maximization = 1)
  int OptDir;
  // optimization direction coefficient (1 for min, -1 for max)
  double OptDirCoeff;
  // true if opt dir has been set manually by the user (don't overwrite
  // with data from problem)
  bool ManualOptDir;

  // current solution
  vector<double> x;
  double f;

  // best solution found
  vector<double> xstar;
  double fstar;

  // starting point
  vector<double> StartingPoint;

  // variable bounds
  vector<double> lb;
  vector<double> ub;

  // is the current best solution feasible?
  bool IsFeasible;

  // optimal objective function value
  double OptimalObjVal;

  // local solver name
  string LocalSolverName;

  // local solver pointer
  Solver* LocalSolver;

  // maximum number of local solutions (corresponds to maximum
  // number of cuts in the local solver)
  int MaxLocalSolutions;

  // maximum running time in seconds (disabled by default)
  double MaxRunningTime;

  // tolerance below which an improvement is not considered so
  double ImprovingTolerance;
  
  // maximum number of non-improving outer iterations (triggers termination)
  int MaxNonImprovingIterations;

  // maximum tabu list length
  int MaxListLength;

  // minimum ball radius
  double MinRadius;

  // maximum ball radius
  double MaxRadius;

  // progressive increase or decrease in radius size
  // (negative DeltaRadius means: start from MaxRadius and decrease tabu
  // region)
  double DeltaRadius;
  
  // whether to print out msgs or not
  bool Quiet;

  // parameters
  Parameters ParameterBlob;

  // name of main solver
  string MainSolver;
  
  // whether rose called from ampl or not
  bool AmplFlag;

};

#endif
