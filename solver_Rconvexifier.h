/*
** Name:       solver_Rconvexifier.h
** Author:     Leo Liberti & Sonia Cafieri
** Source:     GNU C++
** Purpose:    solver Rconvexifier header
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    080606 work started
*/

#ifndef _ROSESOLVERRconvexifierH
#define _ROSESOLVERRconvexifierH

#include "solver.h"

#define ROSEINFINITY 1e30

class RconvexifierSolver : public virtual Solver {

 public:
  RconvexifierSolver();
  ~RconvexifierSolver();

  // data access
  void SetProblem(Problem* p); 
  string GetName(void) { return TheName; }
  Problem* GetProblem(void) { return InProb; }
  Problem* GetProblem(double *dlb, double *dub, map<int, int> var_map);
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

  Expression addLbConstVar(double coeff, int index, vector <double> *vlb, vector <double> *vub);
  Expression addUbConstVar(double coeff, int index, vector <double> *vlb, vector <double> *vub);
  //map to constant variables
  map<int, int> lbconstmap;
  map<int, int> ubconstmap;
  int nconvar; //number of constant variables


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
  int NumberOfObjectives;

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

  map<int,int> VarIDToLocalIndex;
  // optimal objective function value
  double OptimalObjVal;

  // parameters
  Parameters ParameterBlob;

  // be quiet
  bool Quiet;

  // name of main solver
  string MainSolver;
  
  // true if called from ampl
  bool AmplFlag;

  int flagconvex;

  // output .ros file name
  string OutFile;

  // output an AMPL file mapping added variable indices to 
  // corresponding relaxation constraint names
  string OutFileMaps;

  // if true, disables simplifier and outputs 
  // Rconvexifier_out.run and Rprintmod_out.mod
  bool OutAmpl;
};

#endif
