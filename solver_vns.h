/*
** Name:       solver_vns.h
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    solver VNS header 
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    050921 work started
*/

#ifndef _ROSESOLVERVNSH
#define _ROSESOLVERVNSH

#include "solver.h"

#define ROSEINFINITY 1e30

class VNSSolver : public virtual Solver {

 public:
  VNSSolver();
  ~VNSSolver();

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

  // sample in k-th neighbourhood
  // neighbourhoodtype = 0: hyperrectangular shells
  // neighbourhoodtype = 1: hyperrectangles
  void GetRandomPointInNeighbourhood(std::vector<double>& x,
				     std::vector<double>& c,
				     std::vector<double>& l,
				     std::vector<double>& u,
				     int k,
				     int neighbourhoodtype);

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
  
  // best solution
  vector<double> xstar;
  double fstar;

  // variable bounds
  vector<double> vlb;
  vector<double> vub;

  // variable integrality
  vector<bool> integrality;

  // optimal objective function value
  double OptimalObjVal;

  // parameters
  Parameters ParameterBlob;

  // maximum running time
  double MaxRunningTime;
  
  // Kmin, Kmax
  int Kmin;
  int Kmax;

  // samples
  int Samples;

  // is problem mixed-integer?
  bool IsMINLP;

  // name of local solver
  string LocalSolverName;
  
  // local continuous NLP solver pointer
  Solver* LocalSolver;

  // name of local MINLP solver
  string MINLPLocalSolverName;

  // local MINLP solver pointer
  Solver* MINLPLocalSolver;

  // type of neighbourhood
  // (0 = hyperrectangular shells, 1 = hyperrectangles (filled))
  int NeighbourhoodType;

  // be quiet
  bool Quiet;

  // to deal with quiet in local solver
  bool savequiet;

  // relax integer variables?
  bool Relaxed;

  // number of added cuts
  int NumberOfCuts;
  
  // maximum number of allowed added cuts (defaults to 20)
  int MaxNumberOfCuts;

  // do we need cuts?
  bool NeedsCuts;

  // name of main solver
  string MainSolver;
  
  // whether rose called from ampl or not
  bool AmplFlag;

  // warm up local solver first?
  bool WarmUp;

};

#endif
